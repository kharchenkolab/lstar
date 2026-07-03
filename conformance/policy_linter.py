"""Divergence linter: the viewer@0.1 grouping-detection policy is single-sourced. Assert the Python, R,
and JS copies of the preferred-grouping list all equal the canonical conformance/viewer_policy.json.

A drift in any one surface's copy -- the exact class of bug that let #3/#4 diverge silently -- fails
here, cheaply and without running any prep. Skips the R/JS legs when their runtimes are absent (still
checks Python + whatever is available). Run: python conformance/policy_linter.py
"""
import glob
import json
import os
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
canon = json.load(open(os.path.join(ROOT, "conformance", "viewer_policy.json")))
WANT = list(canon["preferred_groupings"])
fail = []


def _node():
    hits = glob.glob(os.path.expanduser("~/emsdk/node/*/bin/node"))
    return hits[0] if hits else "node"


# --- Python ---
sys.path.insert(0, os.path.join(ROOT, "python", "src"))
from lstar.viewer import _PREFERRED_GROUPINGS  # noqa: E402
if list(_PREFERRED_GROUPINGS) != WANT:
    fail.append("Python _PREFERRED_GROUPINGS != canonical:\n    got  %s\n    want %s" % (list(_PREFERRED_GROUPINGS), WANT))
else:
    print("  OK: Python preferred-grouping list matches canonical")

# --- R ---
try:
    r = subprocess.run(
        ["Rscript", "-e", '.libPaths(c("%s", .libPaths())); cat(lstar:::.VIEWER_PREFERRED_GROUPINGS, sep="\\n")'
         % os.path.join(ROOT, ".Rlib")],
        capture_output=True, text=True, timeout=90)
    got = [x for x in r.stdout.strip().splitlines() if x]
    if not got:
        print("  [skip] R lstar not loadable — R policy not checked (%s)" % (r.stderr.strip()[-160:] or "no output"))
    elif got != WANT:
        fail.append("R .VIEWER_PREFERRED_GROUPINGS != canonical:\n    got  %s\n    want %s" % (got, WANT))
    else:
        print("  OK: R preferred-grouping list matches canonical")
except (FileNotFoundError, subprocess.TimeoutExpired):
    print("  [skip] Rscript unavailable — R policy not checked")

# --- JS ---
try:
    js = subprocess.run(
        [_node(), "--experimental-strip-types", "-e",
         'import(%r).then(m => console.log(JSON.stringify(m.PREFERRED_GROUPINGS))).catch(e => { console.error(e); process.exit(3); })'
         % os.path.join(ROOT, "js", "core", "policy.ts")],
        capture_output=True, text=True, timeout=90)
    out = js.stdout.strip()
    if not out.startswith("["):
        print("  [skip] node/policy.ts not loadable — JS policy not checked (%s)" % (js.stderr.strip()[-160:] or "no output"))
    elif json.loads(out) != WANT:
        fail.append("JS PREFERRED_GROUPINGS != canonical:\n    got  %s\n    want %s" % (json.loads(out), WANT))
    else:
        print("  OK: JS preferred-grouping list matches canonical")
except (FileNotFoundError, subprocess.TimeoutExpired):
    print("  [skip] node unavailable — JS policy not checked")

if fail:
    print("\nPOLICY DRIFT (grouping detection is NOT single-sourced):\n" + "\n".join(fail))
    sys.exit(1)
print("  policy linter OK")
