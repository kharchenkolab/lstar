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
