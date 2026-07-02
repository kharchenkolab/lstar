// @lstar/core — read L* Zarr stores in the browser/Node and run the viewer queries.
export { openLstar, LstarDataset } from "./reader.ts";
export type { LstarStore, AxisMeta, FieldMeta } from "./reader.ts";
export { LstarView, Crossfilter, scalarToRGBA } from "./view.ts";
export type { ColStats, Metadata } from "./view.ts";
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
