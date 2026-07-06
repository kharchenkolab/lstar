#!/usr/bin/env bash
# ============================================================================================
# retest_local.sh -- the LOCAL full-corpus retest orchestrator (run rarely; see RETEST.md).
#
# One entry point for the recurring loop: sweep the profiles across the *large local corpus*,
# check that the committed SYNTHETIC CI fixtures still faithfully mirror the real data, and emit
# a TRIAGE telling you exactly what to update (a profile, a synthetic, or a doc). It does NOT edit
# synthetics or docs itself -- those need judgment; it tells you precisely where.
#
# The local corpus is gitignored (large/many) and lives in `testdata/` + the `.Rlib`/Bioconductor/
# SeuratData caches. Every step is TOLERANT: a missing dataset, an absent package, or a single
# crashing dataset degrades to SKIP/LOADERR, never aborts the run. So this is safe to run as-is;
# it covers whatever you have cached.
#
# Usage:
#   bash conformance/sweep/retest_local.sh                 # guard + all cached-data sweeps + triage
#   RETEST_HEAVY=1   bash conformance/sweep/retest_local.sh # also run the GB-scale auto-installing
#                                                           # SeuratData sweeps (seurat/integration/refs)
#   RETEST_GATE=1    bash conformance/sweep/retest_local.sh # also run the full conformance gate first
#                                                           # (run.sh, REAL loaders -- the unit/round-trip suite)
#   RETEST_ONLY="mudata spatial"  bash ...                 # restrict to named sweeps (space-separated keys)
# Acquiring the datasets: see README.md ("How to acquire the local datasets"). The recurring
# directive (when to run, how to triage, how to update synth + docs): see RETEST.md.
# ============================================================================================
set -uo pipefail
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$ROOT"
RLIB="$ROOT/.Rlib"
PY="PYTHONPATH=$ROOT/python/src python3"
HEAVY="${RETEST_HEAVY:-0}"; GATE="${RETEST_GATE:-0}"; ONLY="${RETEST_ONLY:-}"
LOG="/tmp/retest_local_$(date +%s 2>/dev/null || echo run).log"
rm -f /tmp/sweep_*.tsv 2>/dev/null || true
hr(){ printf '%.0s-' {1..92}; echo; }
want(){ [ -z "$ONLY" ] && return 0; case " $ONLY " in *" $1 "*) return 0;; *) return 1;; esac; }
run(){ # run() <label> <command...>  -- announce, run tolerantly, record rc
  local label="$1"; shift
  echo "  ▸ $label"; ( eval "$*" ) >>"$LOG" 2>&1 && echo "    ok" || echo "    (exit $? -- see $LOG; sweeps record SKIP/LOADERR inside their .tsv)"
}

echo "== lstar LOCAL corpus retest =="
echo "   root=$ROOT  heavy=$HEAVY gate=$GATE only='${ONLY:-all cached}'  log=$LOG"
hr

# --- (0) optional: the full conformance gate on REAL loaders (unit + round-trip suite) -------------
if [ "$GATE" = "1" ]; then
  echo "[0] conformance gate (run.sh, REAL loaders -- LSTAR_SYNTHETIC_CORPUS unset)"
  ( unset LSTAR_SYNTHETIC_CORPUS; bash conformance/run.sh ) 2>&1 | tee -a "$LOG" | tail -6
  hr
fi

# --- (1) faithfulness guard: do the committed SYNTHETIC fixtures still cover the REAL structure? ----
# This is the synth-update trigger. A gap here ("synthetic missing <role/axis/subtype>") means a real
# structure is no longer represented in synth.py -> update synth.py, NOT a profile bug.
echo "[1] synthetic-faithfulness guard (real structure must be a subset of the synthetic fixture)"
GUARD=PASS
if want guard; then
  if eval "$PY python/tests/test_synth_faithful.py" 2>&1 | tee -a "$LOG" | sed 's/^/    /'; then :; else GUARD=GAP; fi
else echo "    (skipped by RETEST_ONLY)"; fi
hr

# --- (2) corpus sweeps (each tolerant: missing data -> the sweep SKIPs / prints "none cached") -------
echo "[2] corpus sweeps -> /tmp/sweep_*.tsv  (FAIL/VALIDATE-ERR = profile bug; LOADERR/SKIP = data/pkg absent)"
# Python sweeps: read cached testdata/ with lstar's own deps; subprocess-isolated internally where noted.
want v3           && run "v3 format on real data (v2<->v3 across C++/Py/JS/R)"  "bash conformance/v3_corpus.sh"
want anndata      && run "anndata (scanpy + local atlases + velocity + CITE)" "$PY conformance/sweep/sweep_anndata.py"
want velocity     && run "velocity (scVelo spliced/unspliced)"               "$PY conformance/sweep/sweep_velocity.py"
want mudata       && run "mudata (.h5mu: minipbcite + 10x CITE/multiome)"    "$PY conformance/sweep/sweep_mudata.py"
want spatial      && run "spatial AnnData (Visium + squidpy imaging)"        "$PY conformance/sweep/sweep_spatial.py"
want perturbation && run "perturbation (scPerturb backed read)"             "$PY conformance/sweep/sweep_perturbation.py"
# R sweeps over cached data (no GB installs): scRNAseq (61 SCEs, isolated), 10x CITE/multiome, conos, spatial.
want scrnaseq     && run "scRNAseq (61 Bioconductor SCEs, isolated)"         "bash conformance/sweep/sweep_scrnaseq_driver.sh"
want citeseq_10x  && run "Seurat 10x CITE-seq (.h5 RNA+ADT)"                 "R_LIBS=$RLIB Rscript conformance/sweep/sweep_citeseq_10x.R"
want multiome_10x && run "Seurat 10x multiome (.h5 RNA+ChromatinAssay)"      "R_LIBS=$RLIB Rscript conformance/sweep/sweep_multiome_10x.R"
want conos        && run "Conos (local .rds collections)"                    "R_LIBS=$RLIB Rscript conformance/sweep/sweep_conos.R"
want spatial_r    && run "Seurat spatial (stxBrain + ssHippo, auto-installs)" "R_LIBS=$RLIB Rscript conformance/sweep/sweep_spatial.R"
# Heavy SeuratData sweeps install GBs on first run -> opt in with RETEST_HEAVY=1.
if [ "$HEAVY" = "1" ]; then
  want seurat       && run "SeuratData enumeration (heavy installs)"          "R_LIBS=$RLIB Rscript conformance/sweep/install_and_sweep_seurat.R"
  want integration  && run "SeuratData integration AS collections"            "R_LIBS=$RLIB Rscript conformance/sweep/sweep_integration.R"
  want seurat_refs  && run "SeuratData Azimuth ref atlases"                   "R_LIBS=$RLIB Rscript conformance/sweep/sweep_seurat_refs.R"
else
  echo "  (heavy SeuratData sweeps skipped; set RETEST_HEAVY=1 to include seurat/integration/refs)"
fi
hr

# --- (3) aggregate + triage --------------------------------------------------------------------------
echo "[3] consolidated results"
total_fail=0; total_validate=0; rows=""
shopt -s nullglob
for tsv in /tmp/sweep_*.tsv; do
  nm=$(basename "$tsv" .tsv | sed 's/^sweep_//')
  # status is a column whose values are PASS/FAIL/VALIDATE-ERR/LOADERR/SKIP/CRASH/TIMEOUT (header line skipped)
  # grep -c always prints a number (0 on no match); no `|| echo 0` -- that would double-print "0"
  p=$(grep -cwE "PASS" "$tsv" 2>/dev/null); p=${p:-0}
  f=$(grep -cwE "FAIL" "$tsv" 2>/dev/null); f=${f:-0}
  v=$(grep -cwE "VALIDATE-ERR|VALIDATE_ERR|VALERR" "$tsv" 2>/dev/null); v=${v:-0}
  o=$(grep -cwE "LOADERR|SKIP|CRASH|TIMEOUT" "$tsv" 2>/dev/null); o=${o:-0}
  total_fail=$((total_fail + f)); total_validate=$((total_validate + v))
  printf "    %-16s PASS=%-4s FAIL=%-3s VALIDATE-ERR=%-3s other(skip/load/crash)=%-3s\n" "$nm" "$p" "$f" "$v" "$o"
  rows="$rows$nm:$f:$v;"
done
[ -z "$rows" ] && echo "    (no /tmp/sweep_*.tsv produced -- no cached datasets? see README.md acquisition)"
hr

echo "== TRIAGE -- the recurring loop (full directive: conformance/sweep/RETEST.md) =="
ACTION=0
if [ "$GUARD" = "GAP" ]; then
  echo "  [synth] FAITHFULNESS GAP -> a real structure is no longer covered by a synthetic fixture."
  echo "          Update python/tests/synth.py to add the missing role/axis/subtype, then re-run [1]."
  echo "          (The guard prints exactly which key drifted, e.g. 'subtypes: synthetic missing [spatial]'.)"
  ACTION=1
fi
if [ "$total_fail" -gt 0 ] || [ "$total_validate" -gt 0 ]; then
  echo "  [profile] $total_fail FAIL + $total_validate VALIDATE-ERR across sweeps -> a PROFILE BUG on real data."
  echo "            Inspect the offending rows: grep -E 'FAIL|VALIDATE' /tmp/sweep_*.tsv"
  echo "            Fix the profile (profiles/*.py | R/R/profile_*.R), add a fixture that reproduces it,"
  echo "            then re-sweep. New real structure that now passes also belongs in a synthetic fixture."
  ACTION=1
fi
echo "  [docs] After any change, refresh the docs that quantify coverage:"
echo "         - conformance/sweep/REPORT.md   (the tally: per-modality PASS/FAIL + each finding)"
echo "         - conformance/sweep/CATALOG.md  (the corpus map: dataset -> what it exercises -> CI synthetic)"
echo "         - SUPPORT.md                    (the user-facing matrix: status / real / CI-synth columns)"
echo "         A finding flips from 'gap' to 'Fixed (<commit>)'; a newly-covered case gets its row + a CI fixture."
[ "$ACTION" = "0" ] && echo "  [clean] guard passed and no profile FAIL/VALIDATE-ERR -- only refresh docs if the corpus changed."
hr
echo "full log: $LOG   (sweeps: /tmp/sweep_*.tsv)"
# exit non-zero on anything that needs action -- a profile regression (FAIL/VALIDATE-ERR) OR a
# faithfulness gap (synthetics drifted) -- so this can gate a pre-release / pre-merge check.
if [ "$total_fail" -gt 0 ] || [ "$total_validate" -gt 0 ] || [ "$GUARD" = "GAP" ]; then exit 1; fi
exit 0
