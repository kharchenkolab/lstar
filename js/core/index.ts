// @lstar/core — read L* Zarr stores in the browser/Node and run the viewer queries.
export { openLstar, LstarDataset } from "./reader.ts";
export type { LstarStore, AxisMeta, FieldMeta } from "./reader.ts";
export { LstarView, Crossfilter, scalarToRGBA } from "./view.ts";
export type { ColStats, Metadata } from "./view.ts";
// viewer-compute recipe (single-sourced): pure reductions over CSC measure + codes + cell-sets, shared by
// extend (prep), LstarView (live), and pagoda3. Callers read the measure (fieldAsCsc / csrRows) and reduce.
export { colStats, overdispersionScore, overdispersionFromStats, groupSufficientStats, groupSizes, markers, deAvsB, kernels } from "./compute.ts";
export type { CscMeasure } from "./compute.ts";
export { writeStore, addToStore } from "./writer.ts";
export type { LstarWritableStore, AxisSpec, FieldSpec, DatasetSpec,
  Compressor, WriteOptions, AuxSpec, AuxArraySpec } from "./writer.ts";
// viewer@0.1 optimization (JS twin of Python extend_for_viewer): precompute the viewer navigators into a store.
export { extendForViewer } from "./extend.ts";
export type { ExtendOptions } from "./extend.ts";
// pick the count basis by state (raw preferred; counts=/basis= override) — shared by extend + the viewer prep.
export { selectCountsBasis } from "./basis.ts";
export type { CountsBasis } from "./basis.ts";
// HttpStore adds byte-range (`Range`) reads on top of a plain fetch, enabling the reader's sub-chunk
// fast path over HTTP/CDN. Browser code can also use zarrita's FetchStore (no range fast path) or any
// store implementing the optional `getRange`; NodeFSStore is in ./node-store.ts (imported separately
// so the browser bundle never pulls in node:fs).
export { HttpStore } from "./http-store.ts";
// Single-file `.lstar.zarr.zip` (STORED): ZipStore reads a chunk by ONE range read into the archive
// (HTTP Range or file pread), no extraction — the reason the archive must be STORED. httpZipSource wraps
// a URL as the byte source; packStoredZip writes an in-memory store to a STORED zip. Node file/dir
// helpers (nodeFileSource, openLstarZip, packStoredZipDir, writeStoreZip) live in ./zip-node.ts.
export { ZipStore, httpZipSource, packStoredZip, crc32, readZipCentralDir } from "./zip.ts";
export type { ByteSource } from "./zip.ts";
