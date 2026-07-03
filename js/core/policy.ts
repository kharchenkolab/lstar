// Canonical viewer@0.1 grouping-detection policy for the JS surface. MUST match Python
// (lstar.viewer._PREFERRED_GROUPINGS) and R (.VIEWER_PREFERRED_GROUPINGS); enforced against
// conformance/viewer_policy.json by conformance/policy_linter.py. Kept in its own module with NO
// runtime deps so the linter can import it without loading the WASM kernels.
export const PREFERRED_GROUPINGS = [
  "leiden", "cluster", "clusters", "cell_type", "celltype", "cell_types",
  "louvain", "seurat_clusters", "annotation", "cluster_label",
];
export const MIN_GROUPS = 2;
export const MAX_GROUPS = 60;

// Preference rank of a label name: index of the first preferred term that is a substring of the
// lowercased name (so matches sort first, by list position); non-matches rank last. Ties are broken by
// the caller (alphabetical). Mirrors Python `_rank` / R `.viewer_grouping_rank`.
export function groupingRank(name: string): number {
  const low = name.toLowerCase();
  for (let i = 0; i < PREFERRED_GROUPINGS.length; i++) if (low.includes(PREFERRED_GROUPINGS[i])) return i;
  return PREFERRED_GROUPINGS.length;
}

// These three mirror viewer_policy.json (lognorm_name_fallback, preferred_embeddings, hilbert_grid) and
// are enforced against it by conformance/policy_linter.py — single-sourced here so basis.ts / extend.ts
// don't hardcode them.
export const LOGNORM_NAMES = ["X", "data", "logcounts"];   // lognorm measure-name fallback
export const PREFERRED_EMBEDDINGS = ["umap"];
export const HILBERT_GRID = 1024;

// Preference rank of an embedding name (mirrors groupingRank): index of the first preferred term that is
// a substring of the lowercased name; non-matches rank last. Ties broken by the caller (alphabetical).
export function embeddingRank(name: string): number {
  const low = name.toLowerCase();
  for (let i = 0; i < PREFERRED_EMBEDDINGS.length; i++) if (low.includes(PREFERRED_EMBEDDINGS[i])) return i;
  return PREFERRED_EMBEDDINGS.length;
}
