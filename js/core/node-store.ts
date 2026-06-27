// A node filesystem store (for tests / local CLI use). In the browser, use zarrita's FetchStore
// (HTTP/S3/CDN), ZipFileStore, or a File System Access adapter instead — the reader only needs the
// `get(key)` contract.
import * as fs from "node:fs/promises";
import * as path from "node:path";

import type { LstarStore } from "./reader.ts";

export class NodeFSStore implements LstarStore {
  root: string;
  constructor(root: string) { this.root = root; }
  async get(key: string): Promise<Uint8Array | undefined> {
    try {
      return new Uint8Array(await fs.readFile(path.join(this.root, key)));
    } catch (e: any) {
      if (e && e.code === "ENOENT") return undefined;
      throw e;
    }
  }
  /** Read bytes [start, end) of one object (a positional read), or undefined if the object is absent.
   * Enables the reader's sub-chunk byte-range fast path without loading the whole array. */
  async getRange(key: string, start: number, end: number): Promise<Uint8Array | undefined> {
    let fh: fs.FileHandle | undefined;
    try {
      fh = await fs.open(path.join(this.root, key), "r");
      const len = Math.max(0, end - start);
      const buf = new Uint8Array(len);
      if (len === 0) return buf;
      const { bytesRead } = await fh.read(buf, 0, len, start);
      return bytesRead === len ? buf : buf.subarray(0, bytesRead);
    } catch (e: any) {
      if (e && e.code === "ENOENT") return undefined;
      throw e;
    } finally {
      await fh?.close();
    }
  }
  /** Write one object (mkdir -p the parent). The write side of the `get` contract used by writer.ts. */
  async set(key: string, value: Uint8Array): Promise<void> {
    const fp = path.join(this.root, key);
    await fs.mkdir(path.dirname(fp), { recursive: true });
    await fs.writeFile(fp, value);
  }
  /** Remove one object if present (used to drop stale consolidated metadata after a write). */
  async delete(key: string): Promise<void> {
    await fs.rm(path.join(this.root, key), { force: true });
  }
}
