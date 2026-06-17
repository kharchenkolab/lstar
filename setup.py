"""Build the optional C++ accelerator (lstar._accel) over the header-only libstar core.

Metadata lives in pyproject.toml; this file only defines the compiled extension and makes it
*optional* -- if the compiler/OpenMP are unavailable, the package still installs and runs on the
pure-Python fallback (lstar.has_accel() will report False).
"""
import os
import sys

from setuptools import Extension, setup
from setuptools.command.build_ext import build_ext


def _libomp_prefix():
    """Locate Homebrew's keg-only libomp (arm64 /opt/homebrew, x86_64 /usr/local) so omp.h / -lomp
    resolve -- without this the macOS accelerator build fails with 'omp.h file not found'."""
    import subprocess
    try:
        p = subprocess.run(["brew", "--prefix", "libomp"], capture_output=True, text=True)
        cand = p.stdout.strip()
        if cand and os.path.isdir(cand):
            return cand
    except Exception:  # noqa: BLE001  (brew may be absent off-CI)
        pass
    for cand in ("/opt/homebrew/opt/libomp", "/usr/local/opt/libomp"):
        if os.path.isdir(cand):
            return cand
    return None


def _ext():
    import pybind11
    include_dirs = [pybind11.get_include(), "core/include"]
    if sys.platform == "darwin":
        cflags = ["-std=c++17", "-O3", "-Xpreprocessor", "-fopenmp"]
        lflags = ["-lomp"]
        omp = _libomp_prefix()                     # keg-only: add its include/lib paths explicitly
        if omp:
            include_dirs.append(os.path.join(omp, "include"))
            cflags.append("-I" + os.path.join(omp, "include"))
            lflags.append("-L" + os.path.join(omp, "lib"))
    elif sys.platform == "win32":
        cflags = ["/std:c++17", "/O2", "/openmp", "/EHsc"]
        lflags = []
    else:
        cflags = ["-std=c++17", "-O3", "-fopenmp"]
        lflags = ["-fopenmp"]
    return Extension(
        "lstar._accel",
        sources=["python/src/lstar/_accel.cpp"],
        include_dirs=include_dirs,
        extra_compile_args=cflags,
        extra_link_args=lflags,
        language="c++",
    )


class OptionalBuildExt(build_ext):
    """Never fail the install because the accelerator didn't compile -- fall back to pure Python."""

    def run(self):
        try:
            super().run()
        except Exception as e:  # noqa: BLE001
            sys.stderr.write("lstar: C++ accelerator not built (%s); using pure-Python\n" % e)

    def build_extension(self, ext):
        try:
            super().build_extension(ext)
        except Exception as e:  # noqa: BLE001
            sys.stderr.write("lstar: skipping accelerator %s (%s)\n" % (ext.name, e))


setup(ext_modules=[_ext()], cmdclass={"build_ext": OptionalBuildExt})
