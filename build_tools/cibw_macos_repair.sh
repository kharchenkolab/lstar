#!/usr/bin/env bash
# macOS wheel repair -- make the bundled libomp dyld-coalescible so a second OpenMP runtime never loads.
#
# THE BUG (arm64 SIGSEGV in extend_for_viewer / any OpenMP kernel):
#   `delocate-wheel` bundles libomp into lstar/.dylibs and rewrites _accel's dependency to the
#   PATH-SPECIFIC install name `@loader_path/.dylibs/libomp.dylib`. dyld deduplicates loaded images by
#   install name; a concrete @loader_path path does NOT match the `@rpath/libomp.dylib` that other
#   scientific packages (scikit-learn, numba, some scipy builds) use -- so when one of those has already
#   loaded its libomp, a SECOND copy maps into the process. Two llvm-openmp runtimes then coexist and their
#   __kmp_* worker/suspend machinery cross-binds between the copies; on macOS arm64 that corrupts thread
#   state and crashes (EXC_BAD_ACCESS). KMP_DUPLICATE_LIB_OK=TRUE does not help; OMP_NUM_THREADS=1 masks it.
#
# THE FIX (option B -- keep the bundle, restore coalescing):
#   Rename _accel's dependency back to `@rpath/libomp.dylib`. Now dyld REUSES a host libomp already loaded
#   under that name (single runtime, no crash), while an `@loader_path/.dylibs` rpath still resolves our
#   bundled copy when no host libomp is present (a bare `pip install lstar-sc` still imports). This mirrors
#   how the ecosystem avoids multi-libomp crashes: everyone names it `@rpath/libomp.dylib` so dyld loads one.
#
# Verify with build_tools/macos_libomp_repro.py (forces a 2nd libomp via scikit-learn, then runs the
# kernels). cibuildwheel invokes this as the macOS `repair-wheel-command` (see pyproject.toml).
set -euo pipefail
wheel="$1"; dest_dir="$2"; archs="${3:-$(uname -m)}"

# 1) normal delocate: bundle libomp into lstar/.dylibs, ad-hoc sign, emit the repaired wheel to $dest_dir.
delocate-wheel --require-archs "$archs" -w "$dest_dir" -v "$wheel"
repaired="$dest_dir/$(basename "$wheel")"

# 2) rewrite the _accel -> libomp dependency name inside the wheel. Unpack/pack via the `wheel` tool so the
#    RECORD hashes are recomputed after we mutate + re-sign the .so (a raw `zip` update would leave a stale
#    RECORD and pip would reject the install).
python3 -m pip install -q --disable-pip-version-check wheel
work="$(mktemp -d)"
python3 -m wheel unpack -d "$work" "$repaired"
so="$(find "$work" -name '_accel*.so' -print -quit)"
[ -n "$so" ] || { echo "cibw_macos_repair: no _accel*.so found in $repaired" >&2; exit 1; }

if otool -L "$so" | grep -q '@loader_path/.dylibs/libomp.dylib'; then
  install_name_tool -change @loader_path/.dylibs/libomp.dylib @rpath/libomp.dylib "$so"
  otool -l "$so" | grep -q 'path @loader_path/.dylibs' \
    || install_name_tool -add_rpath @loader_path/.dylibs "$so"
  codesign --remove-signature "$so" 2>/dev/null || true
  codesign --force --sign - "$so"                 # arm64 requires a valid signature; install_name_tool voids it
  rm -f "$repaired"
  ( cd "$work"/* && python3 -m wheel pack -d "$dest_dir" . )   # repack -> recomputes RECORD, restores filename
  echo "cibw_macos_repair: renamed _accel libomp dependency -> @rpath/libomp.dylib (dyld-coalescible)"
  otool -L "$(find "$dest_dir" -name '*.whl' -newer "$work" 2>/dev/null | head -1)" 2>/dev/null || true
else
  echo "cibw_macos_repair: _accel has no @loader_path libomp dependency; left wheel unchanged"
fi
rm -rf "$work"
