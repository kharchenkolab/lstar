"""Compute-engine selection: fast C++ by default, transparent pure-Python fallback.

The compiled `lstar._accel` extension (a binding over the libstar OpenMP kernels) is used whenever
it imported successfully; otherwise the pure-Python implementation runs. Users choose nothing -- a
`pip install` of a wheel gives the fast path automatically, and an environment without the built
extension still works, just slower. `engine=` on the public calls and the `LSTAR_ENGINE` env var
are escape hatches for benchmarking/debugging, not part of the normal path.
"""
import os

try:
    from . import _accel  # compiled extension (present in wheels / after a source build)
    _HAVE_ACCEL = True
except Exception:          # noqa: BLE001 -- any import failure means "no accelerator", fall back
    _accel = None
    _HAVE_ACCEL = False


def has_accel():
    """True if the compiled C++ accelerator is available (the fast path is in use by default)."""
    return _HAVE_ACCEL


def resolve_engine(engine):
    """Resolve engine ('auto'|None|'c++'|'python') to 'c++' or 'python'.

    'auto'/None picks C++ when available (honoring LSTAR_ENGINE for an override); 'c++' requires the
    extension and raises if absent; 'python' forces the reference path.
    """
    if engine in (None, "auto"):
        env = os.environ.get("LSTAR_ENGINE", "").strip().lower()
        if env in ("python", "py"):
            return "python"
        if env in ("c++", "cpp", "cxx"):
            engine = "c++"
        else:
            return "c++" if _HAVE_ACCEL else "python"
    if engine in ("c++", "cpp", "cxx"):
        if not _HAVE_ACCEL:
            raise RuntimeError(
                "engine='c++' was requested but the lstar._accel extension is not available. "
                "Install a wheel that bundles the compiled core, or use engine='python'.")
        return "c++"
    if engine in ("python", "py"):
        return "python"
    raise ValueError("unknown engine %r (use 'auto', 'c++', or 'python')" % (engine,))


def show_config():
    """Print which compute path lstar will use -- the first thing to check if it seems slow."""
    if _HAVE_ACCEL:
        omp = getattr(_accel, "openmp", False)
        print("lstar compute: C++ accelerator (libstar) ACTIVE")
        print("  OpenMP: %s%s" % (omp, "  (max threads %d)" % _accel.max_threads() if omp else ""))
        print("  default engine: 'c++'  (override with engine=... or LSTAR_ENGINE=python)")
    else:
        print("lstar compute: pure-Python fallback (the _accel extension did not import)")
        print("  install a binary wheel for the C++ accelerator; default engine: 'python'")
