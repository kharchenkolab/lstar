#!/usr/bin/env bash
# Install a LOW-deployment-target libomp for the macOS wheel build.
#
# Homebrew's libomp is built for the *runner's* macOS (13/14), so bundling it into a wheel forces a recent
# macOS minimum -- delocate then rejects the wheel ("Library dependencies do not satisfy target ..."). The
# conda-forge `llvm-openmp` package is built for an OLD target (macOS 10.x/11.0), so a wheel that bundles it
# runs on macOS >= MACOSX_DEPLOYMENT_TARGET (we set 11.0). Mirrors the scikit-learn / scipy approach.
#
# Extracts to $LSTAR_LIBOMP (default ~/libomp): <prefix>/include/omp.h, <prefix>/lib/libomp.dylib. setup.py
# picks it up via LSTAR_LIBOMP. No conda/brew runtime needed -- just a download + untar.
set -euo pipefail

case "$(uname -m)" in
  arm64)  SUBDIR=osx-arm64 ;;
  x86_64) SUBDIR=osx-64 ;;
  *) echo "cibw_macos_openmp: unsupported arch $(uname -m)" >&2; exit 1 ;;
esac
DEST="${LSTAR_LIBOMP:-$HOME/libomp}"

python3 - "$SUBDIR" "$DEST" <<'PY'
import io, json, os, sys, tarfile, urllib.request
subdir, dest = sys.argv[1], sys.argv[2]
# Resolve the newest conda-forge llvm-openmp .tar.bz2 for this arch via the anaconda.org API (no hash to
# guess, no conda needed). .tar.bz2 is a plain tar (unlike the newer .conda zstd zip), so stdlib untars it.
meta = json.load(urllib.request.urlopen("https://api.anaconda.org/package/conda-forge/llvm-openmp"))
cands = [f for f in meta["files"]
         if f["attrs"].get("subdir") == subdir and f["basename"].endswith(".tar.bz2")]
if not cands:
    sys.exit("cibw_macos_openmp: no .tar.bz2 llvm-openmp build for %s" % subdir)
cands.sort(key=lambda f: f["attrs"].get("timestamp", 0))
f = cands[-1]
url = f["download_url"]
url = ("https:" + url) if url.startswith("//") else url
sys.stderr.write("cibw_macos_openmp: %s\n" % url)
buf = io.BytesIO(urllib.request.urlopen(url).read())
os.makedirs(dest, exist_ok=True)
with tarfile.open(fileobj=buf, mode="r:bz2") as t:
    t.extractall(dest)                       # -> dest/include/omp.h, dest/lib/libomp.dylib
assert os.path.isfile(os.path.join(dest, "lib", "libomp.dylib")), "libomp.dylib missing after extract"
assert os.path.isfile(os.path.join(dest, "include", "omp.h")), "omp.h missing after extract"
sys.stderr.write("cibw_macos_openmp: installed low-target libomp -> %s\n" % dest)
PY
