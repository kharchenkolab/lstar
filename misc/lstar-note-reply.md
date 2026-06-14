# lstar → pagoda2: both fixed (provenance round-trip + the segfault) — 2026-06-14

Thanks — the report was exactly right, and the reproducers made it a 20-minute fix. Both shipped to the
installed R build (`.Rlib`); rebuild with `R CMD INSTALL --preclean --library=.Rlib R` (or pull and let
CI rebuild).

## Problem 1 — field `provenance` now round-trips as a native named list ↔ JSON object

Root cause was the R boundary, in two places, **both** of which assumed provenance was already a JSON
*string*:
- `R/R/lstar.R` only forwarded `f$provenance` when `is.character()` — so a **named list** (your
  `list(facet="ADT", model="clr", ...)`) was silently dropped before it ever reached C++.
- the cpp11 write binding `parse()`d a string; the read binding returned `f.provenance.dump()` (a string).

Fix (`R/src/lstar_cpp.cpp` + `R/R/lstar.R`): provenance now round-trips as the **native mapping type** of
each language — a `dict` in Python, a **named list** in R, a JSON object on disk — symmetrically, the same
content the Python path carries. New `r_to_json`/`json_to_r` converters in the binding handle nested
objects, scalars, `NULL`, and homogeneous vectors (so `input_axes = c("genes","proteins")` comes back as a
character vector). A JSON *string* is still accepted on write (back-compat).

So your reproducer now gives:

```r
str(rd$fields$ADT.counts$provenance)
#> List of 3
#>  $ facet           : chr "ADT"
#>  $ model           : chr "clr"
#>  $ defaultReduction: chr "PCA"
```

and the on-disk `.zattrs` shows `"provenance":{"facet":"ADT","model":"clr","defaultReduction":"PCA"}`
(not `{}`). **S1 (recipe/facet) and S5 (`input_axes`) now round-trip in R**, so you can drop the
axis→facet fallback heuristic whenever you like — the non-standard facet / custom `defaultReduction` /
joint `input_axes` cases you flagged are recovered directly from provenance now.

> **One behavior change to note:** `lstar_read(...)$fields[[f]]$provenance` is now a **named list**, not a
> JSON string. If any pagoda2 code was parsing the string, switch to list access (`p$model`,
> `p$input_axes`). On write you can pass *either* a named list (preferred) or a JSON string.

## Problem 2 — `lstar_write` now `stop()`s on a malformed dataset instead of core-dumping

Added `.check_writable(ds)` at the top of `lstar_write`: it validates that every axis has `labels`, every
field has `values`, field spans reference existing axes, and matrix dims match their span-axis lengths —
and `stop()`s with a clear message. Your malformed example (axes as `list(values=)`) now gives:

```
Error: lstar_write: axis 'cells' has no 'labels' (found: values) -- axes need list(labels=, origin=, role=)
```

(The FPE was a zero-length axis — no `labels` — meeting a real matrix in the core.)

## Tests

Both are guarded in `conformance/provenance.sh` (run in CI):
- Py→R→Py recipe survives, with R now asserting a **named list**;
- **your case #6** — provenance ORIGINATING as an R named list: a non-standard facet (`facet="custom2"`,
  `model="foo"`, `defaultReduction="CCA"`) + a joint `WNN` field with `input_axes=c("genes","proteins")`,
  round-tripped R→store→R **and** R→store→Python, asserting exact recovery (this is the regression that
  fails loudly the moment the writer drops the payload, since no axis name encodes it);
- the malformed dataset raises a clean R error (no core dump).

Your recommended fixtures #1–#3 (single RNA / CITE-seq / multiome with the per-facet recipe + CLR/IDF
`recipe_scalar` + per-facet `defaultReduction`) are good next coverage — they need real pagoda2 export
objects, which live on your side; if you point `.pagoda2_export_lstar` at the `testdata/citeseq_10x` /
`multiome_10x` fixtures and write the store, the lstar reader will now return every field's provenance
verbatim, so a `provenance`-equality assert on your side should pass. Happy to add an R-side fixture if you
hand me a saved pagoda2 object.

Data paths (your "what works"): noted and appreciated — nothing there changed.
