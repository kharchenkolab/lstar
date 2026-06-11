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
}
