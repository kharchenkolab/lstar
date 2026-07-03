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

  // explicit measure — honour it; log1p unless it (or basis=) says the values are already log-normalized.
  if (counts != null) {
    const f = ds.field(counts);
    if (!f) throw new Error(`extendForViewer: counts="${counts}" is not a measure (present cells×genes measures: ${present})`);
    return { field: counts, log1p: basis !== "lognorm" && f.state !== "lognorm" };
  }
  // explicit lognorm basis — a log-normalized measure (by state, else the usual names), used as-is.
  if (basis === "lognorm") {
    const pick = twod.find((n) => ds.field(n)!.state === "lognorm") ?? twod.find((n) => LOGNORM_NAMES.includes(n));
    if (!pick) throw new Error(`extendForViewer: basis="lognorm" but no log-normalized measure found (present: ${present})`);
    return { field: pick, log1p: false };
  }
  // default: a measure named "counts", else any raw-state measure; log1p'd.
  const pick = (twod.includes("counts") ? "counts" : undefined) ?? twod.find((n) => ds.field(n)!.state === "raw");
  if (pick) return { field: pick, log1p: true };
  throw new Error(
    `extendForViewer: not viewer-optimizable — no raw counts measure found (present cells×genes measures: ${present}). ` +
    `Viewer prep needs raw counts; pass counts=<field>, provide a raw-counts measure, or pass ` +
    `basis="lognorm" to prep (approximately) from a log-normalized measure.`);
}
