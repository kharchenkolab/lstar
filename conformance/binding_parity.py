"""Binding-parity tripwire: a compute kernel bound in one surface's binding must be bound in ALL of them,
unless it is explicitly declared surface-scoped below. This catches the class where a core kernel is added
to (or registered in) one binding but not the others -- e.g. `viewer_cell_order` was defined in _accel.cpp
but never registered, so Python silently ran the numpy fallback instead of the OpenMP core; and a kernel
bound only in Python/WASM would slip past an R reader.

Parses the REGISTERED kernel names from each binding (Python `m.def`, WASM `function(...)`, the generated
R cpp11 registration table), normalizes to snake_case, and asserts symmetry modulo the ALLOW list. No
compiler or runtime needed -- pure text. Run: python conformance/binding_parity.py
"""
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _read(*p):
    return open(os.path.join(ROOT, *p)).read()


def _camel_to_snake(s):
    return re.sub(r"(?<!^)(?=[A-Z])", "_", s).lower()


py = set(re.findall(r'm\.def\("([a-z_0-9]+)"', _read("python", "src", "lstar", "_accel.cpp")))
wasm = {_camel_to_snake(n) for n in re.findall(r'function\("([A-Za-z0-9]+)"', _read("js", "wasm", "lstar_wasm.cpp"))}
r = set(re.findall(r'"_lstar_lstar_cpp_([a-z_0-9]+)"', _read("R", "src", "cpp11.cpp")))

# Legitimately surface-scoped bindings (NOT a parity bug). Adding a kernel to only some surfaces requires a
# line here with the reason -- that's the tripwire: an *undeclared* asymmetry fails. See docs/parity.md.
ALLOW = {
    "max_threads": "Python-only thread-count utility",
    "version": "WASM-only build/version string",
    "gzip_compress": "WASM-only codec (Py uses numcodecs; R via C++ core)",
    "csc_to_csr": "encoding flip -- Py/WASM bind it; R uses the Matrix package",
    "csr_to_csc": "encoding flip -- WASM-only (no scipy/Matrix in the browser); Py/R normalize natively",
    "col_mean_var": "Py/WASM bind it directly; R reaches it via stream_col_stats",
    "read": "R-only store IO (Py uses zarr-python, JS uses zarrita)",
    "write": "R-only store IO (Py uses zarr-python, JS uses zarrita)",
    "read_csc_block": "R-only bounded block read (Py/JS use their own zarr layer)",
    "read_csc_cols": "R-only bounded column read (Py/JS use their own zarr layer)",
    "stream_col_stats": "R/C++-only depth-normalized streaming reducer (pagoda2 host)",
    "stream_col_sum_by_group": "R/C++-only streamed pseudobulk (pagoda2 host)",
}

allf = py | wasm | r
fail = []
for k in sorted(allf):
    if k in ALLOW:
        continue
    where = {"Python": k in py, "R": k in r, "WASM": k in wasm}
    if not all(where.values()):
        present = [s for s, v in where.items() if v]
        missing = [s for s, v in where.items() if not v]
        fail.append("kernel '%s' bound in %s but MISSING from %s -- bind it there, or add to ALLOW with a "
                    "reason if it's genuinely surface-scoped" % (k, present, missing))

print("  registered: Python=%d R=%d WASM=%d; %d surface-scoped (allowed)" % (len(py), len(r), len(wasm), len(ALLOW)))
if fail:
    print("\nBINDING-PARITY DRIFT (a shared kernel is bound unevenly):\n" + "\n".join("  " + f for f in fail))
    sys.exit(1)
print("  OK: shared compute kernels bound on all three surfaces: " + ", ".join(sorted(allf - set(ALLOW))))
