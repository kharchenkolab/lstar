# Suggestions for lstar from the pagoda2.1 multimodal work

*Design pressure surfaced while drafting the pagoda2.1 multimodal model (the `facet` redesign +
storage-backed Conos collections). Source: `pagoda2/misc/multimodal_proposal2.md` §7, §8.6.4. Written
2026-06-13. These are **candidate** changes for the lstar side to weigh — none is assumed.*

## Framing — keep lstar's algorithmic scope bounded

The guiding constraint (from the pagoda2 side, but it protects lstar): **lstar is persistence +
interchange + a thin, general kernel layer — not a compute/method engine.** Basic, modality-agnostic
primitives belong in the core; anything recipe- or method-specific belongs in the *user* package
(pagoda2, cacoa, …), with a clean way to keep it there. Every item below is checked against that line.
The one place lstar has already crossed it — "DE / markers / pseudobulk bundles" (`SUPPORT.md:62`) — is
called out in S4.

Status legend mirrors `SUPPORT.md`: `✓` done · `◐` partial · `✗` not yet.

---

## S1 — Persist normalization *recipe parameters* as provenance; do **not** compute (✓ likely already)

**Problem.** pagoda2.1 stores *only* raw counts plus a small normalization recipe (`plain`/`clr`/`tfidf`
+ params); the normalized "analysis view" is derived on the fly and is **never** a stored matrix. For the
interchange to stay equally lean, a pagoda2 → lstar → pagoda2 round-trip must not be forced to either (a)
write a redundant dense `lognorm` measure, or (b) silently lose *which* normalization was intended.

**Non-ask (retracted).** An earlier draft asked lstar to implement `recipe` **virtual fields** it
materializes on read. That was wrong: a computed field isn't persisted data, and making lstar evaluate
normalization is exactly the scope creep S1 wants to avoid.

**Ask (narrow).** Just preserve the recipe *parameters* as a typed metadata record on the **raw**
measure — `depthScale`, log base, winsor caps, `model`, and (for byte-exact reproduction) the per-axis
precomputed scalars (CLR per-row divisor, per-column IDF). lstar's existing `state` (`raw`/`lognorm`/…)
+ `provenance` already accommodate this; the round-trip stays lean, and any kernel (pagoda2's `misc2.cpp`,
`libstar`, or the WASM viewer) recomputes values on demand at a boundary that needs them (AnnData `X` on
export; the viewer). Likely **already supported** via provenance passthrough — please confirm it survives
a round-trip rather than being dropped.

**S1′ (optional convenience).** Agree a *typed* recipe vocabulary (`plain`/`clr`/`tfidf` + params) so the
in-browser WASM kernel (which already runs view kernels, `SUPPORT.md:66`) can apply CLR/TF-IDF from the
params without pagoda2 in the loop. Optional — opaque provenance round-trips faithfully; typing only buys
cross-kernel reuse.

---

## S2 — Collection-level streaming grouped reducer (extend the *primitive*, not a bundle) (✗ / ◐)

**Problem.** The storage-backed Conos path (pagoda2 §8.6) must compute pseudobulk / per-cluster summaries
over an atlas-scale collection without materializing the joint matrix. The natural operation: aggregate
per-sample `counts.{s}` measures into one `(joint-cluster × gene)` result over the **derived union**
`cells` axis, streamed, in bounded memory.

**Ask.** lstar lists a grouped-sum kernel and DE/pseudobulk bundles (`SUPPORT.md:60,62`), but presumably
within a *single* measure. Extend the **grouped-sum primitive** so it walks a collection's per-sample
measures and accumulates into the induced factor axis (rule-2 group axis), with float64-over-float32
accumulation. Frame this as a primitive extension — *not* a new "markers bundle" (see S4).

---

## S3 — Confirm `index` partial coverage composes with a *derived union* axis (◐ — needs a case)

**Problem.** A pagoda2.1 facet membership mask maps 1:1 to an lstar partial-coverage `index`
(`SUPPORT.md:55`, `python/tests/test_partial.py`). At collection scale the index must point into the
**derived union** `cells` axis of a `kind="collection"`, not only an observed axis.

**Ask.** Add a conformance case: a partial measure with `index_axis` = a derived union axis, round-tripped
Py↔C++↔R. If it already works, document it; if not, it's a small fix on the validate/IO path.

---

## S4 — Determinism contract (take it) + optional bounded shared core (decide carefully)

**S4a — the determinism contract (free, no scope; ✓ wants anyway).** A *property*, not new
functionality: reducers are thread-count-invariant and bit-identical — float64 accumulation,
column-parallel, **no** cross-thread reduction (pagoda2 §6.2). lstar already advertises reproducible,
thread-count-controllable reductions (`principles.md §5`); the only ask is to make it a **stated, tested
contract** so a facet's summaries are identical whichever library computes them.

**S4b — optional convergence on one core, with a hard boundary (decision, not a default).** pagoda2's
fused-reducer/zarr-gather work overlaps `libstar`'s `col mean/var · fused depth+log1p · grouped sum`
(`SUPPORT.md:60`), so sharing one tuned, deterministic core is tempting. It is only safe with a boundary
that keeps lstar bounded:

- **In the core (short, closed list):** lazy block IO with deterministic decode order, plus ~4 reduction
  skeletons taking `(raw block, per-row scalar r[], per-col scalar c[], transform enum t)`, where
  `t ∈ {identity, ÷rowscalar, ×colscalar, log1p, winsorize}`. The §6.2 property — a *pure per-entry
  function of precomputed scalars* — is simultaneously what makes the kernel thread-invariant and what
  keeps method-specificity out.
- **In the user package (everything else):** the code that *computes* the scalars (CLR divisor, IDF,
  winsor caps, variance model) lives in pagoda2; it hands the core a vector + enum, and the core never
  learns what "CLR" is.
- **Two-tier extension, neither tier grows lstar:** (1) expressible as per-entry(scalars, enum) → call the
  existing skeleton, no new core code; (2) anything else (e.g. a cross-entry transform) → use `libstar`'s
  block reader + **header-only template** instantiation and write the loop in the user package.
- **Decision rule:** converge only while pagoda2's kernel set fits the scalars+enum boundary (today it
  does); the first kernel that *can't* be expressed that way is the signal to keep it in `misc2.cpp`, not
  to widen `libstar`.

**Caution (the line lstar has already crossed).** "DE / markers / pseudobulk bundles" (`SUPPORT.md:62`)
are past the primitive line. The method-specific parts (which test, which statistic, specificity metrics)
are candidates to push *out* of the core toward the user package or a profile, on the same principle that
bounds S4b.

---

## S5 — Provenance for named-product input facets (✗, minor)

**Problem.** A joint reduction (WNN/MOFA/totalVI) is a *named product* over the shared cells axis with the
contributing facets as provenance — pagoda2 stores `reductions[["WNN"]]` with `facets=c("RNA","ADT")`.

**Ask.** lstar provenance should carry an explicit input-feature-axis / facet list on an embedding so the
named product's facet set round-trips losslessly. lstar already round-trips the *shape* of these (WNN
weights as cell measures + joint graph; MOFA shared-factor scores + per-mod loadings, `SUPPORT.md:107-109`)
— this just records *which facets* fed it.

---

## Priority for the lstar side

| item | what it is | effort | recommendation |
|---|---|---|---|
| **S4a** | determinism contract (stated + tested) | tiny | take it — lstar wants it regardless |
| **S1** | recipe params as provenance (confirm round-trip) | tiny | confirm; likely already works |
| **S3** | `index` over a derived union axis (conformance case) | small | add the case |
| **S5** | facet-set provenance on named products | small | additive |
| **S2** | collection-level grouped-sum across per-sample measures | medium | the storage-backed Conos path needs it |
| **S4b** | converge pagoda2/lstar on one C++ core | large | **a real architectural decision** — only with the S4b boundary; otherwise keep `misc2.cpp` and take only S4a |

The architecturally significant choice is **S4b** (one shared deterministic kernel core, bounded by the
scalars+enum line). Everything else is either free (S4a, S1) or additive (S2/S3/S5). The retracted item is
the old "virtual `recipe` field" — lstar should *not* compute; it should only persist the recipe params.
