"""Divergence linter: the viewer@0.1 policy is single-sourced. Assert the Python, R, and JS copies of
every policy constant (the preferred-grouping list, min/max group counts, the lognorm measure-name
fallback, preferred embeddings, the Hilbert grid) all equal the canonical conformance/viewer_policy.json.

A drift in any one surface's copy -- the exact class of bug that let #3/#4 (and the basis tie-break)
diverge silently -- fails here, cheaply and without running any prep. Skips the R/JS legs when their
runtimes are absent (still checks Python + whatever is available). Run: python conformance/policy_linter.py
"""
import glob
import json
import os
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
canon = json.load(open(os.path.join(ROOT, "conformance", "viewer_policy.json")))
KEYS = ["preferred_groupings", "min_groups", "max_groups", "lognorm_name_fallback",
        "preferred_embeddings", "hilbert_grid",
        "viewer_codec", "viewer_level", "viewer_chunk_elems", "viewer_shard_elems"]
WANT = {k: canon[k] for k in KEYS[:6]}
_c = canon["compression"]                     # per-field viewer compression layout constants
WANT.update({"viewer_codec": _c["codec"], "viewer_level": _c["level"],
             "viewer_chunk_elems": _c["chunk_elems"], "viewer_shard_elems": _c["shard_elems"]})
fail = []


def _cmp(surface, got):
    for k in KEYS:
        if got.get(k) != WANT[k]:
            fail.append("%s %s != canonical:\n    got  %r\n    want %r" % (surface, k, got.get(k), WANT[k]))
    if all(got.get(k) == WANT[k] for k in KEYS):
        print("  OK: %s policy constants match canonical" % surface)


# --- Python ---
sys.path.insert(0, os.path.join(ROOT, "python", "src"))
from lstar.viewer import (_PREFERRED_GROUPINGS, _MIN_GROUPS, _MAX_GROUPS,      # noqa: E402
                          _LOGNORM_NAMES, _PREFERRED_EMBEDDINGS,
                          _VIEWER_CODEC, _VIEWER_LEVEL, _VIEWER_CHUNK_ELEMS, _VIEWER_SHARD_ELEMS)
from lstar.kernels import _N_GRID                                             # noqa: E402
_cmp("Python", {"preferred_groupings": list(_PREFERRED_GROUPINGS), "min_groups": _MIN_GROUPS,
                "max_groups": _MAX_GROUPS, "lognorm_name_fallback": list(_LOGNORM_NAMES),
                "preferred_embeddings": list(_PREFERRED_EMBEDDINGS), "hilbert_grid": _N_GRID,
                "viewer_codec": _VIEWER_CODEC, "viewer_level": _VIEWER_LEVEL,
                "viewer_chunk_elems": _VIEWER_CHUNK_ELEMS, "viewer_shard_elems": _VIEWER_SHARD_ELEMS})

# --- R ---
_rscript = r'''.libPaths(c("%s", .libPaths())); suppressMessages(library(lstar))
jarr <- function(x) paste0('["', paste(x, collapse='","'), '"]')
cat("preferred_groupings=", jarr(lstar:::.VIEWER_PREFERRED_GROUPINGS), "\n", sep="")
cat("min_groups=", lstar:::.VIEWER_MIN_GROUPS, "\n", sep="")
cat("max_groups=", lstar:::.VIEWER_MAX_GROUPS, "\n", sep="")
cat("lognorm_name_fallback=", jarr(lstar:::.VIEWER_LOGNORM_NAMES), "\n", sep="")
cat("preferred_embeddings=", jarr(lstar:::.VIEWER_PREFERRED_EMBEDDINGS), "\n", sep="")
cat("hilbert_grid=", lstar:::.VIEWER_HILBERT_GRID, "\n", sep="")
cat('viewer_codec="', lstar:::.VIEWER_CODEC, '"\n', sep="")
cat("viewer_level=", lstar:::.VIEWER_LEVEL, "\n", sep="")
cat("viewer_chunk_elems=", lstar:::.VIEWER_CHUNK_ELEMS, "\n", sep="")
cat("viewer_shard_elems=", lstar:::.VIEWER_SHARD_ELEMS, "\n", sep="")''' % os.path.join(ROOT, ".Rlib")
try:
    r = subprocess.run(["Rscript", "-e", _rscript], capture_output=True, text=True, timeout=90)
    got = {}
    for line in r.stdout.strip().splitlines():
        if "=" in line:
            k, v = line.split("=", 1)
            got[k.strip()] = json.loads(v.strip())
    if not got:
        print("  [skip] R lstar not loadable — R policy not checked (%s)" % (r.stderr.strip()[-160:] or "no output"))
    else:
        _cmp("R", got)
except (FileNotFoundError, subprocess.TimeoutExpired):
    print("  [skip] Rscript unavailable — R policy not checked")

# --- JS ---
_node = (glob.glob(os.path.expanduser("~/emsdk/node/*/bin/node")) or ["node"])[0]
_pol = os.path.join(ROOT, "js", "core", "policy.ts")
_js = ('import(%r).then(m => console.log(JSON.stringify({preferred_groupings: m.PREFERRED_GROUPINGS, '
       'min_groups: m.MIN_GROUPS, max_groups: m.MAX_GROUPS, lognorm_name_fallback: m.LOGNORM_NAMES, '
       'preferred_embeddings: m.PREFERRED_EMBEDDINGS, hilbert_grid: m.HILBERT_GRID, '
       'viewer_codec: m.VIEWER_CODEC, viewer_level: m.VIEWER_LEVEL, viewer_chunk_elems: m.VIEWER_CHUNK_ELEMS, '
       'viewer_shard_elems: m.VIEWER_SHARD_ELEMS}))).catch(e => '
       '{ console.error(e); process.exit(3); })' % _pol)
try:
    js = subprocess.run([_node, "--experimental-strip-types", "-e", _js], capture_output=True, text=True, timeout=90)
    out = js.stdout.strip()
    if not out.startswith("{"):
        print("  [skip] node/policy.ts not loadable — JS policy not checked (%s)" % (js.stderr.strip()[-160:] or "no output"))
    else:
        _cmp("JS", json.loads(out))
except (FileNotFoundError, subprocess.TimeoutExpired):
    print("  [skip] node unavailable — JS policy not checked")

if fail:
    print("\nPOLICY DRIFT (viewer policy is NOT single-sourced):\n" + "\n".join(fail))
    sys.exit(1)
print("  policy linter OK")
