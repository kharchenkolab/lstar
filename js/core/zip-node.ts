// Node-only helpers for single-file `.lstar.zarr.zip` stores: a file `ByteSource` (so a local zip is
// seek-into-readable), a directory packer, and read/write conveniences. Kept out of ./zip.ts so the
// browser bundle never pulls in node:fs; browser code uses ZipStore + httpZipSource / packStoredZip.
import * as fsp from "node:fs/promises";
import * as os from "node:os";
import * as path from "node:path";

import { openLstar, type LstarDataset } from "./reader.ts";
import { NodeFSStore } from "./node-store.ts";
import { writeStore, type DatasetSpec, type WriteOptions } from "./writer.ts";
import { ZipStore, packStoredZip, readZipCentralDir } from "./zip.ts";
import type { ByteSource } from "./zip.ts";

/** A `ByteSource` over a local file: size via `stat`, ranges via positional reads (`pread`). */
export function nodeFileSource(filePath: string): ByteSource {
  return {
    async size() { return (await fsp.stat(filePath)).size; },
    async range(start, end) {
      const fh = await fsp.open(filePath, "r");
      try {
        const len = Math.max(0, end - start);
        const buf = new Uint8Array(len);
        if (len === 0) return buf;
        const { bytesRead } = await fh.read(buf, 0, len, start);
        return bytesRead === len ? buf : buf.subarray(0, bytesRead);
      } finally { await fh.close(); }
    },
  };
}

/** Open a local `.lstar.zarr.zip` as an L* dataset (seek-into-zip, no extraction). */
export async function openLstarZip(zipPath: string): Promise<LstarDataset> {
  return openLstar(await ZipStore.open(nodeFileSource(zipPath), zipPath));
}

/** Pack a directory store into ONE STORED file (ZIP64-aware). Mirrors the Python/C++ packers. */
export async function packStoredZipDir(dir: string, zipPath: string): Promise<void> {
  const entries: Array<[string, Uint8Array]> = [];
  async function walk(rel: string): Promise<void> {
    for (const de of await fsp.readdir(path.join(dir, rel), { withFileTypes: true })) {
      const childRel = rel ? rel + "/" + de.name : de.name;
      if (de.isDirectory()) await walk(childRel);
      else entries.push([childRel, new Uint8Array(await fsp.readFile(path.join(dir, childRel)))]);
    }
  }
  await walk("");
  await fsp.writeFile(zipPath, packStoredZip(entries));
}

/** Extract every STORED entry of a `.lstar.zarr.zip` into `dir` (a copy per entry, no decompression) —
 * e.g. to obtain a writable directory store from a zip before extending it. Rejects DEFLATE (via
 * readZipCentralDir's guard). */
export async function extractZipToDir(zipPath: string, dir: string): Promise<void> {
  const src = nodeFileSource(zipPath);
  const idx = await readZipCentralDir(src, zipPath);
  const zs = await ZipStore.open(src, zipPath);
  for (const key of idx.keys()) {
    if (!key || key.endsWith("/")) continue;
    const bytes = await zs.get(key);
    const fp = path.join(dir, key);
    await fsp.mkdir(path.dirname(fp), { recursive: true });
    await fsp.writeFile(fp, bytes!);
  }
}

/** Write a complete L* store as ONE STORED `.lstar.zarr.zip` (writes a temp dir, then packs it). */
export async function writeStoreZip(zipPath: string, ds: DatasetSpec, opts?: WriteOptions): Promise<void> {
  const tmp = await fsp.mkdtemp(path.join(os.tmpdir(), "lstar-zip-"));
  try {
    await writeStore(new NodeFSStore(tmp), ds, opts);
    await packStoredZipDir(tmp, zipPath);
  } finally {
    await fsp.rm(tmp, { recursive: true, force: true });
  }
}
