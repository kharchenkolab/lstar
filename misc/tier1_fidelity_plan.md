# Tier 1 fidelity — implementation plan (P0)

*The two load-bearing Tier-1 gaps from `format_coverage.md`, to land before P3 (DE). Both are
**baseline correctness** (silent round-trip corruption today), and both must round-trip across **all
four** implementations — Python, C++, R, JS — with conformance gates. Sequenced B → A; each phase is
gated on a cross-language round-trip before moving on.*

## Part B — nullable / extension dtypes (values + validity mask)

**Problem.** pandas nullable `Int64`/`boolean`/`string` columns (values + a null mask) are coerced to
float-NaN (losing integer-ness) or to `"<NA>"` strings. Distinct from categorical `-1`-missing (which
P1 fixed) and from float NaN (which already encodes missing). The fix is an explicit **validity mask**.

**Design.** A field may carry an optional **`mask`** companion (`uint8`, `1 = missing`) beside its
`values`, for the `dense` (integer/bool) and `utf8` (string) encodings. Float fields keep using NaN (no
mask). The mask is the one mechanism; it does not change the value encoding.

- **on disk:** `fields/<name>/mask` (uint8, same length as values) + `lstar.nullable=true`. Absent mask
  ⇒ no nulls (back-compatible: existing stores read unchanged).
- **model.py:** `Field.mask` (optional ndarray); `add_field(..., mask=)`; a `Masked` value wrapper is
  *not* needed — keep `values` raw + a sibling `mask`. `_infer_*` unchanged.
- **zarr_io.py:** write `mask` when present; read it back onto `Field.mask`; lazy path materializes it.
- **validate.py:** `mask` length == axis length; dtype uint8/bool.
- **C++ (`lstar.hpp`):** `Field` gains `NdArray mask; bool has_mask`; read/write the `mask` array +
  `nullable` flag. Re-vendor to R.
- **R (`lstar_cpp.cpp`, `lstar.R`):** read → reconstruct an NA-bearing vector (integer `NA`, logical
  `NA`, character `NA`); write → split a vector with `NA` into values (fill) + mask. Marshal `mask`.
- **JS (`reader.ts`/`view.ts`):** `FieldMeta.nullable`; `fieldDense`/`fieldStrings` expose a `mask`;
  `view.metadata()` carries it so a null renders distinctly (not 0).
- **anndata profile:** `_by_dtype_series` detects `pd.api.types.is_extension_array_dtype` (Int/boolean/
  string) → emit raw values + mask; write-back rebuilds `pd.array(values, dtype=..., mask=...)`.
- **gate:** `conformance/nullable.sh` — a masked integer + masked string written in Python read
  back in C++/R (NA reconstructed) and re-written, read in Python with mask byte-identical; +
  `test_nullable.py`; + JS reader assertion in `js.sh`.

## Part A — lossless `uns` / `@misc` passthrough (the `aux/` subtree)

**Problem.** `uns` / `@misc` / `@commands` are recorded **name-only** in `dropped` — the content
(params, colors, variance ratios, dendrograms, DE tables) is lost. Upgrade `dropped` from "record the
name" to **verbatim preserve + reproduce on write**. Foundational: makes round-trips safe *before* we
type anything, and lets recognized structures be promoted out of the tail incrementally.

**Design — a generic, self-describing passthrough subtree** (not an opaque blob: it stays inspectable
and promotable, the doc's stated goal; reuses L*'s own primitives — JSON attrs + typed arrays). A
reserved **`aux/`** group, listed in root meta as `aux:[...]`. Each entry `aux/<ns>` (namespaced, e.g.
`anndata.uns`) holds:

```
aux/<ns>/.zattrs  lstar = { kind:"aux", tree:<JSON>, arrays:[ {id, enc}, … ] }
aux/<ns>/<id>            a typed array (enc="dense")  OR  utf8 bytes+offsets (enc="utf8")
```

`tree` is the nested structure with array leaves replaced by references; every dense/utf8 leaf is listed
flat in `arrays`. **Leaf grammar in `tree`:** JSON scalar (inline) · `{"$array":id}` (→ a dense array)
· `{"$strings":id}` (→ utf8) · `{"$record":{fields:{name:leaf}, length:n}}` (numpy structured array,
e.g. `rank_genes_groups`) · nested object/array (recurse).

This split means **C++/R round-trip is purely mechanical**: copy `tree` JSON verbatim + iterate the
`arrays` manifest reading/writing each leaf — *no tree interpretation*. Only the originating profile
(Python anndata) walks `tree` to reconstruct the live object. So:

- **core serializer (`aux.py`, language-neutral algorithm):** `to_store(obj)->(tree, arrays)` /
  `from_store(tree, arrays)->obj` over: scalar · str · list · dict · ndarray (num/str/bool) · structured
  ndarray. Genuinely unrepresentable leaves (sparse, opaque objects) → still recorded in `dropped`
  (best-effort lossless, honest about the remainder).
- **model.py:** `Dataset.aux: dict[str, Any]` (reconstructed objects, Python-side).
- **zarr_io.py:** write/read the `aux/` group; root meta gains `aux`.
- **C++ (`lstar.hpp`):** `struct Aux { string ns; json tree; vector<ArrayLeaf> arrays; }`;
  `Dataset.aux`; read (tree from `.zattrs`, leaves per manifest) + write. Re-vendor.
- **R:** round-trip `ds$aux[[ns]]` = list(tree=<parsed>, arrays=<named list>); reconstruct-to-R-list is a
  nice-to-have after the round-trip gate.
- **JS (`reader.ts`):** read `aux` → `{tree, arrays}` resolved into a JS object (read-only, for the
  viewer / promotion).
- **anndata profile:** `read_anndata` → `ds.aux["anndata.uns"] = adata.uns` (minus the parts already
  typed: neighbors graphs already in obsp, etc. — keep it simple first: stash whole `uns`, let typed
  fields win on write); `write_anndata` → reconstruct `adata.uns`. Replaces the name-only `dropped`.
- **gate:** `conformance/aux.sh` — (1) a **synthetic** nested tree (dict/list/scalars + num/str/record
  arrays) round-trips Py↔C++↔R byte-identical (core mechanism, format-agnostic); (2) an AnnData with a
  rich `uns` (params, `*_colors`, `pca.variance_ratio`, `dendrogram` linkage, a `rank_genes_groups`
  structured array) round-trips `uns` **exactly** through L* (Py→L*→Py) and survives a C++/R re-write; +
  `test_aux.py`; + JS read assertion.

## Sequencing & gates

```
B nullable (values+mask) ──gate(Py↔C++↔R↔JS round-trip)──▶ A passthrough core (synthetic tree)
   ──gate(Py↔C++↔R byte-identical)──▶ A anndata uns↔aux ──gate(rich-uns exact round-trip)──▶ (P3)
```

- **Out of scope here:** Seurat `@misc`/`@commands` passthrough (the `aux/` core enables it; an R-profile
  follow-on), color-palette *typing*/binding (passthrough already carries `*_colors`; binding to the
  factor axis is a later promotion), `uns` param *typing* (`log1p.base` etc. — promotion, post-P3).
- **Why B first:** smaller, self-contained, establishes the four-language test rhythm; A is the bigger
  foundational piece and reuses that rhythm.
