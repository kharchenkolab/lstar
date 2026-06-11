"""Build the optional C++ accelerator (lstar._accel) over the header-only libstar core.

Metadata lives in pyproject.toml; this file only defines the compiled extension and makes it
*optional* -- if the compiler/OpenMP are unavailable, the package still installs and runs on the
pure-Python fallback (lstar.has_accel() will report False).
"""
import sys

from setuptools import Extension, setup
from setuptools.command.build_ext import build_ext


def _ext():
    import pybind11
    if sys.platform == "darwin":
        cflags = ["-std=c++17", "-O3", "-Xpreprocessor", "-fopenmp"]
        lflags = ["-lomp"]
    elif sys.platform == "win32":
        cflags = ["/std:c++17", "/O2", "/openmp", "/EHsc"]
        lflags = []
    else:
        cflags = ["-std=c++17", "-O3", "-fopenmp"]
        lflags = ["-fopenmp"]
    return Extension(
        "lstar._accel",
        sources=["python/src/lstar/_accel.cpp"],
        include_dirs=[pybind11.get_include(), "core/include"],
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
