"""Build the lstar._accel C++ extension (pybind11 over libstar's core kernels).

Metadata lives in pyproject.toml; this file exists only to declare the compiled extension. The
core header is shared from ../core/include (one kernel, every runtime). OpenMP is enabled when
Homebrew's libomp is present (macOS clang needs -Xpreprocessor -fopenmp + an explicit -lomp);
without it the build is serial but identical in result.

Build in place for development:  python setup.py build_ext --inplace
"""
import os
from pybind11.setup_helpers import Pybind11Extension, build_ext
from setuptools import setup

HERE = os.path.dirname(os.path.abspath(__file__))
CORE_INCLUDE = os.path.join(HERE, "..", "core", "include")

extra_compile = ["-O3"]
extra_link = []
include_dirs = [CORE_INCLUDE]

# OpenMP via Homebrew libomp (optional; the kernels run serial without it).
for omp in ("/opt/homebrew/opt/libomp", "/usr/local/opt/libomp"):
    if os.path.isdir(omp):
        extra_compile += ["-Xpreprocessor", "-fopenmp", "-I%s/include" % omp]
        extra_link += ["-L%s/lib" % omp, "-lomp"]
        break

ext_modules = [
    Pybind11Extension(
        "lstar._accel",
        ["src/lstar/_accel.cpp"],
        include_dirs=include_dirs,
        extra_compile_args=extra_compile,
        extra_link_args=extra_link,
        cxx_std=17,
    )
]

setup(ext_modules=ext_modules, cmdclass={"build_ext": build_ext})
