// @lstar/core — read L* Zarr stores in the browser/Node and run the viewer queries.
export { openLstar, LstarDataset } from "./reader.ts";
export type { LstarStore, AxisMeta, FieldMeta } from "./reader.ts";
export { LstarView, Crossfilter, scalarToRGBA } from "./view.ts";
export type { ColStats, Metadata } from "./view.ts";
export { writeStore, addToStore } from "./writer.ts";
export type { LstarWritableStore, AxisSpec, FieldSpec, DatasetSpec } from "./writer.ts";
// Browser stores come from zarrita directly (FetchStore, etc.); NodeFSStore is in ./node-store.ts
// (imported separately so the browser bundle never pulls in node:fs).
