// Choose the count measure the viewer navigators are built from — by CONTENT/STATE, not the literal
// field name "counts". The JS twin of Python `lstar.viewer._select_counts_basis`, so R / Python / JS
// prep agree on what counts as a raw basis (a converter that named its raw matrix "X" or a modality is
// still viewer-preppable) and fail with the SAME clear message when there is no raw counts — removing
// the whole "converter didn't use the magic name `counts`" failure class across the stack.
//
// Returns { field, log1p }: `log1p` is true for a RAW basis (the kernels log1p it) and false for a
// LOGNORM basis (values used as-is → stats are var-of-lognorm, not var-of-log1p(counts)).
import type { FieldMeta } from "./reader.ts";
import { LOGNORM_NAMES } from "./policy.ts";

export interface CountsBasis { field: string; log1p: boolean; }

// Just the slice of the dataset this needs — so it's trivially unit-testable with a mock.
interface DatasetLike { fieldNames(): string[]; field(name: string): FieldMeta | undefined; }

export function selectCountsBasis(ds: DatasetLike, opts: { counts?: string; basis?: string } = {}): CountsBasis {
  const { counts, basis } = opts;
  // cells×genes measures (+ their state), for selection and for the error messages.
  const twod = ds.fieldNames().filter((n) => {
    const f = ds.field(n);
    return !!f && f.role === "measure" && Array.isArray(f.span) && f.span.length === 2 && String(f.span[0]).startsWith("cells");
  });
  const present = twod.map((n) => `${n}[${ds.field(n)!.state ?? "?"}]`).join(", ") || "(none)";
  const rawPick = () => (twod.includes("counts") ? "counts" : twod.find((n) => ds.field(n)!.state === "raw"));
  // name fallback EXCLUDES a scaled/z-scored measure (a scaled `X` is not lognorm).
  const lognormPick = () => twod.find((n) => ds.field(n)!.state === "lognorm")
    ?? twod.find((n) => LOGNORM_NAMES.includes(n) && ds.field(n)!.state !== "scaled");

  // explicit measure — honour it; log1p unless the values are already log-normalized.
  if (counts != null) {
    const f = ds.field(counts);
    if (!f) throw new Error(`extendForViewer: counts="${counts}" is not a measure (present cells×genes measures: ${present})`);
    return { field: counts, log1p: f.state !== "lognorm" };
  }
  const b = basis ?? "auto";
  if (b === "raw") {
    const pick = rawPick();
    if (!pick) throw new Error(`extendForViewer: basis="raw" but no raw counts measure found (present cells×genes measures: ${present}).`);
    return { field: pick, log1p: true };
  }
  if (b === "lognorm") {
    const pick = lognormPick();
    if (!pick) throw new Error(`extendForViewer: basis="lognorm" but no log-normalized measure found (present cells×genes measures: ${present}).`);
    return { field: pick, log1p: false };
  }
  if (b === "auto") {
    const raw = rawPick();
    if (raw) return { field: raw, log1p: true };            // prefer raw (log1p)
    const ln = lognormPick();
    if (ln) return { field: ln, log1p: false };             // fall back to lognorm (as-is)
    throw new Error(
      `extendForViewer: no raw or log-normalized measure found (present cells×genes measures: ${present}). ` +
      `Viewer prep needs raw counts or log-normalized values; pass counts=<field> to force a measure. ` +
      `(A scaled/z-scored measure cannot be used as a basis.)`);
  }
  throw new Error(`extendForViewer: basis must be "auto", "raw", or "lognorm" (or pass counts=<field>); got "${basis}"`);
}
