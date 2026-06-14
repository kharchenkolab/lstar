"""``lstar convert`` — one-command conversion between single-cell formats, bridged by the L* store.

``lstar convert SRC DST`` detects each format from its path, reads SRC into the L* model, and writes DST
from it. The L* dataset is the universal intermediate, so a conversion is just ``write_Y(read_X(obj))``
and what a target cannot hold is recorded in ``ds.dropped`` (visible, never silently lost).

Each format has two backends: a **native** one (uses the domain package — anndata / SeuratObject / …) and,
where implemented, a **direct** one (lstar's own package-free codec over a base engine: h5py for HDF5,
base R for ``.rds``). ``--backend auto`` (default) prefers native when the package is importable and falls
back to the direct codec otherwise; when neither can handle something, lstar says exactly what to install.

Format detection (override with ``--from`` / ``--to``):
  ``.h5ad`` → anndata · ``.h5mu`` → mudata · ``.rds`` → rds (Seurat/SCE, sniffed R-side) ·
  ``.lstar.zarr`` / ``.zarr`` / a Zarr directory → store
"""
from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
import tempfile


class ConvertError(Exception):
    """A user-facing conversion error (bad path, undetectable format, unsupported route)."""


class NeedsPackage(ConvertError):
    """A wall the package-free path can't cross — names the package + install command to handle it.

    Raised at dispatch (no backend can handle the format) or inside a direct backend (it hit something
    only the native package can represent, e.g. an external-pointer-backed matrix)."""

    def __init__(self, thing: str, package: str, install_cmd: str):
        self.thing, self.package, self.install_cmd = thing, package, install_cmd
        super().__init__(
            f"{thing} needs the '{package}' package. Install it and re-run — lstar uses it "
            f"automatically when present:\n    {install_cmd}")


# Which language materializes each format.
_PY = {"anndata", "mudata", "store"}
_R = {"seurat", "sce", "rds"}

# Native backend per format: (python_pkg | None, r_pkg | None, install_cmd | None).
_PKG = {
    "anndata": ("anndata", None, "pip install anndata"),
    "mudata":  ("mudata", None, "pip install mudata"),
    "seurat":  (None, "SeuratObject", "Rscript -e 'install.packages(\"SeuratObject\")'"),
    "sce":     (None, "SingleCellExperiment", "Rscript -e 'BiocManager::install(\"SingleCellExperiment\")'"),
    "rds":     (None, "SeuratObject", "Rscript -e 'install.packages(\"SeuratObject\")'"),
    "store":   (None, None, None),
}

# Package-free (direct) backends — registered as the Tier-A loops land them; empty == "native only".
_DIRECT_PY_READ: dict = {}     # fmt -> (src) -> Dataset
_DIRECT_PY_WRITE: dict = {}    # fmt -> (ds, dst) -> None
_DIRECT_R_READ: dict = {}      # fmt -> (src, bridge_store) -> None  (writes the store)
_DIRECT_R_WRITE: dict = {}     # fmt -> (bridge_store, dst) -> None


def _direct_anndata_read(src):                      # lazy: only imports h5py/the codec when used
    from .profiles.anndata_direct import read_h5ad_direct
    return read_h5ad_direct(src)


def _direct_anndata_write(ds, dst):
    from .profiles.anndata_direct import write_h5ad_direct
    write_h5ad_direct(ds, dst)


_DIRECT_PY_READ["anndata"] = _direct_anndata_read
_DIRECT_PY_WRITE["anndata"] = _direct_anndata_write

# Extension → format (longest/most-specific first).
_EXT = [(".h5ad", "anndata"), (".h5mu", "mudata"), (".rds", "rds"),
        (".lstar.zarr", "store"), (".zarr", "store")]
# --from/--to aliases.
_ALIAS = {"h5ad": "anndata", "ad": "anndata", "adata": "anndata",
          "h5mu": "mudata", "md": "mudata",
          "zarr": "store", "lstar": "store", "lstar.zarr": "store",
          "rds": "rds", "seurat": "seurat", "sce": "sce",
          "singlecellexperiment": "sce"}


def detect_format(path: str, explicit: str | None = None) -> str:
    """Resolve the L* format name for *path*, honoring an explicit ``--from``/``--to`` override."""
    if explicit:
        f = explicit.lower()
        return _ALIAS.get(f, f)
    p = path.lower()
    for ext, fmt in _EXT:
        if p.endswith(ext):
            return fmt
    if os.path.isdir(path) and any(
            os.path.exists(os.path.join(path, m)) for m in (".zgroup", ".zattrs", ".zmetadata")):
        return "store"
    raise ConvertError(
        f"cannot detect the format of {path!r} from its name — pass --from/--to "
        f"(anndata | mudata | store | seurat | sce)")


# ── backend dispatch ────────────────────────────────────────────────────────────────────────────────

def _has_direct(fmt: str, direction: str) -> bool:
    if fmt in _R:
        return fmt in (_DIRECT_R_READ if direction == "read" else _DIRECT_R_WRITE)
    return fmt in (_DIRECT_PY_READ if direction == "read" else _DIRECT_PY_WRITE)


def _native_available(fmt: str) -> bool:
    pkg_py, pkg_r, _ = _PKG[fmt]
    if pkg_py is not None:
        import importlib.util
        return importlib.util.find_spec(pkg_py) is not None
    if pkg_r is not None:
        if fmt == "rds":                              # a bare .rds source: Seurat *or* SCE class suffices
            return _r_pkg_available("SeuratObject") or _r_pkg_available("SingleCellExperiment")
        return _r_pkg_available(pkg_r)
    return True                                       # store: always native


def _choose_backend(fmt: str, pref: str, direction: str) -> str:
    """Return 'native' or 'direct' for *fmt* under *pref* (auto|native|direct); raise with guidance."""
    if fmt == "store":
        return "native"
    has_direct = _has_direct(fmt, direction)
    pkg_py, pkg_r, install = _PKG[fmt]
    if pref == "direct":
        if not has_direct:
            raise ConvertError(f"no package-free backend for {fmt!r} yet — use --backend native (or auto)")
        return "direct"
    native_ok = _native_available(fmt)
    if pref == "native":
        if native_ok:
            return "native"
        raise NeedsPackage(f"the native {fmt} backend", pkg_py or pkg_r, install)
    # auto: prefer native, fall back to the package-free codec, else say exactly what to install
    if native_ok:
        return "native"
    if has_direct:
        return "direct"
    raise NeedsPackage(f"reading/writing {fmt}", pkg_py or pkg_r, install)


# ── the R engine seam (system Rscript today; webR / R-in-WASM can plug in here for JS/Node hosts) ──────

def run_r_driver(src: str, *args) -> "subprocess.CompletedProcess":
    """Run an R driver (source string) under the configured engine and return the completed process.
    This is the seam where webR can later drop in as an alternate engine; today it is system ``Rscript``."""
    rscript = os.environ.get("LSTAR_RSCRIPT", "Rscript")
    fh = tempfile.NamedTemporaryFile("w", suffix=".R", delete=False)
    fh.write(src)
    fh.close()
    try:
        return subprocess.run([rscript, fh.name, *map(str, args)],
                              stdin=subprocess.DEVNULL, capture_output=True, text=True)
    except FileNotFoundError:
        raise ConvertError(
            "Rscript not found — Seurat/SCE conversion needs R with the lstar package "
            "(set LSTAR_RSCRIPT / LSTAR_RLIB if they're not on the default path)")
    finally:
        os.unlink(fh.name)


_R_AVAIL_CACHE: dict = {}


def _r_pkg_available(pkg: str) -> bool:
    if not pkg:
        return True
    if pkg in _R_AVAIL_CACHE:
        return _R_AVAIL_CACHE[pkg]
    src = ('rlib <- Sys.getenv("LSTAR_RLIB", ""); if (nzchar(rlib)) .libPaths(c(rlib, .libPaths()));'
           f' cat(if (requireNamespace("{pkg}", quietly = TRUE)) "YES" else "NO")')
    try:
        proc = run_r_driver(src)
        ok = "YES" in (proc.stdout or "")
    except ConvertError:
        ok = False                                    # no Rscript at all
    _R_AVAIL_CACHE[pkg] = ok
    return ok


# An R driver run via `Rscript <file>` (a temp file, not `-e` — R's `-e` buffer is ~8 KB and silently
# truncates). It bridges Seurat/SCE <-> the L* store via the native class packages. LSTAR_RLIB prepends
# the lstar R library path.
_R_DRIVER = r'''
args  <- commandArgs(trailingOnly = TRUE)
mode  <- args[1]; path <- args[2]; store <- args[3]
fmt   <- if (length(args) >= 4) args[4] else ""
rlib  <- Sys.getenv("LSTAR_RLIB", "")
if (nzchar(rlib)) .libPaths(c(rlib, .libPaths()))
suppressMessages({
  library(lstar)
  if (requireNamespace("SeuratObject", quietly = TRUE)) library(SeuratObject)
  if (requireNamespace("SingleCellExperiment", quietly = TRUE)) library(SingleCellExperiment)
})
if (identical(mode, "to_store")) {                 # R source (.rds) -> L* store
  obj <- readRDS(path)
  ds  <- if (inherits(obj, "SingleCellExperiment")) read_sce(obj) else read_seurat(obj)
  lstar_write(ds, store)
} else {                                           # L* store -> R target (.rds)
  ds  <- lstar_read(store)
  obj <- if (identical(fmt, "sce")) write_sce(ds) else write_seurat(ds)
  saveRDS(obj, path)
}
cat("LSTAR_R_OK\n")
'''


def _run_r(mode: str, path: str, store: str, fmt: str = "") -> None:
    """Run the native R bridge driver for one leg; raise a clean error on any R failure."""
    proc = run_r_driver(_R_DRIVER, mode, path, store, fmt)
    if proc.returncode != 0 or "LSTAR_R_OK" not in proc.stdout:
        tail = "\n".join((proc.stderr or proc.stdout).strip().splitlines()[-8:])
        raise ConvertError(f"R conversion step failed ({mode}):\n{tail}")


# Package-free (direct) Seurat reader: base R only — readRDS (which auto-loads the class package if it is
# installed, else leaves the class undefined with slots reachable via attr()) + lstar's `.read_seurat_direct`
# (S4 slot-walking, no SeuratObject accessors) + lstar_write. lstar itself only Imports Matrix.
_R_DRIVER_DIRECT_READ = r'''
args <- commandArgs(trailingOnly = TRUE); path <- args[2]; store <- args[3]
rlib <- Sys.getenv("LSTAR_RLIB", ""); if (nzchar(rlib)) .libPaths(c(rlib, .libPaths()))
suppressMessages(library(lstar))
obj <- readRDS(path)
ds  <- lstar:::.read_seurat_direct(obj)
lstar_write(ds, store)
cat("LSTAR_R_OK\n")
'''


def _direct_seurat_read(src: str, bridge: str) -> None:
    proc = run_r_driver(_R_DRIVER_DIRECT_READ, "seurat", src, bridge)
    if proc.returncode != 0 or "LSTAR_R_OK" not in proc.stdout:
        tail = "\n".join((proc.stderr or proc.stdout).strip().splitlines()[-8:])
        raise ConvertError(f"package-free Seurat read failed (--backend direct):\n{tail}")


_DIRECT_R_READ["seurat"] = _direct_seurat_read       # .rds Seurat read with base R only (no SeuratObject)
_DIRECT_R_READ["rds"] = _direct_seurat_read


# ── Python-side load / emit (native or direct, per dispatch) ──────────────────────────────────────────

def _load_py(src: str, fmt: str, backend: str = "auto"):
    """Read a Python-side source into an L* :class:`Dataset`."""
    if not os.path.exists(src):
        raise ConvertError(f"source not found: {src}")
    import lstar
    if fmt == "store":
        return lstar.read(src)
    if _choose_backend(fmt, backend, "read") == "direct":
        return _DIRECT_PY_READ[fmt](src)
    if fmt == "anndata":
        import anndata as ad
        from lstar.profiles.anndata import read_anndata
        return read_anndata(ad.read_h5ad(src))
    if fmt == "mudata":
        import mudata as md
        from lstar.profiles.mudata import read_mudata
        return read_mudata(md.read_h5mu(src))
    raise ConvertError(f"{fmt!r} is not a Python-side source format")


def _emit_py(ds, dst: str, fmt: str, backend: str = "auto") -> None:
    """Write an L* :class:`Dataset` out to a Python-side target."""
    import lstar
    if fmt == "store":
        lstar.write(ds, dst)
        return
    if _choose_backend(fmt, backend, "write") == "direct":
        _DIRECT_PY_WRITE[fmt](ds, dst)
        return
    if fmt == "anndata":
        from lstar.profiles.anndata import write_anndata
        write_anndata(ds).write_h5ad(dst)
        return
    if fmt == "mudata":
        from lstar.profiles.mudata import write_mudata
        write_mudata(ds).write(dst)
        return
    raise ConvertError(f"{fmt!r} is not a Python-side target format")


# ── R-side legs (native bridge or direct base-R codec, per dispatch) ──────────────────────────────────

def _r_read_to_store(src: str, fmt: str, bridge: str, backend: str) -> None:
    if _choose_backend(fmt, backend, "read") == "direct":
        _DIRECT_R_READ[fmt](src, bridge)
    else:
        _run_r("to_store", src, bridge)


def _r_write_from_store(bridge: str, dst: str, fmt: str, backend: str) -> None:
    if _choose_backend(fmt, backend, "write") == "direct":
        _DIRECT_R_WRITE[fmt](bridge, dst)
    else:
        _run_r("from_store", dst, bridge, "sce" if fmt == "sce" else "seurat")


def _bridge_store() -> str:
    return os.path.join(tempfile.mkdtemp(prefix="lstar-convert-"), "bridge.lstar.zarr")


def _read_dataset(src: str, ff: str, backend: str = "auto"):
    """Read any source into an L* :class:`Dataset` — Python in-process, or R (Seurat/SCE) via a temp
    bridge store. Used by ``inspect`` (read-only); ``convert`` keeps its own write-coupled bridging."""
    if ff not in _R:
        return _load_py(src, ff, backend)
    if not os.path.exists(src):
        raise ConvertError(f"source not found: {src}")
    import lstar
    bridge = _bridge_store()
    try:
        _r_read_to_store(src, ff, bridge, backend)
        return lstar.read(bridge)                    # arrays land in memory; the bridge can be removed
    finally:
        shutil.rmtree(os.path.dirname(bridge), ignore_errors=True)


def convert(src: str, dst: str, from_fmt: str | None = None, to_fmt: str | None = None,
            backend: str = "auto"):
    """Convert *src* → *dst*. Returns ``(ds, from_fmt, to_fmt)`` with the bridging L* dataset.

    Python-side ↔ Python-side runs in-process; any leg in R (Seurat/SCE) is bridged through a temporary
    on-disk L* store. *backend* (auto|native|direct) selects the native package vs. lstar's package-free
    codec per leg."""
    import lstar
    ff = detect_format(src, from_fmt)
    tf = detect_format(dst, to_fmt)
    if tf == "rds":
        tf = "seurat"                              # a bare .rds *target* defaults to Seurat (use --to sce)
    src_r, dst_r = ff in _R, tf in _R
    if not src_r and not dst_r:                    # in-process Python path
        ds = _load_py(src, ff, backend)
        _emit_py(ds, dst, tf, backend)
        return ds, ff, tf
    if not os.path.exists(src):
        raise ConvertError(f"source not found: {src}")
    bridge = _bridge_store()                        # cross-language: bridge through a temp L* store
    try:
        if src_r:
            _r_read_to_store(src, ff, bridge, backend)
            ds = lstar.read(bridge)
        else:
            ds = _load_py(src, ff, backend)
            _emit_py(ds, bridge, "store", backend)
        if dst_r:
            _r_write_from_store(bridge, dst, tf, backend)
        else:
            _emit_py(ds, dst, tf, backend)
        return ds, ff, tf
    finally:
        shutil.rmtree(os.path.dirname(bridge), ignore_errors=True)


def _print_summary(ds, src: str, ff: str, dst: str | None, tf: str | None) -> None:
    """A short human summary of what crossed (the default; --report prints the full fidelity report)."""
    head = f"lstar: {os.path.basename(src)} ({ff})"
    head += f"  →  {os.path.basename(dst)} ({tf})" if dst else "  (inspect)"
    print(head)
    print(f"  kind: {ds.kind}   profiles: {', '.join(ds.profiles) or '-'}")
    print(f"  axes: {len(ds.axes)} ({', '.join(list(ds.axes)[:6])}"
          f"{', …' if len(ds.axes) > 6 else ''})")
    print(f"  fields: {len(ds.fields)}")
    dl = list(ds.dropped)
    if dl:
        shown = "; ".join(str(x) for x in dl[:4]) + (" …" if len(dl) > 4 else "")
        print(f"  dropped (not representable in L*): {len(dl)} — {shown}")
    else:
        print("  dropped (not representable in L*): none")


def build_report(ds, src: str, ff: str, dst: str | None = None, tf: str | None = None) -> dict:
    """A structured fidelity report: what crossed (axes + fields with role/state/span/provenance) and
    what could not be represented (``dropped``). Language-neutral — rebuildable from any L* store."""
    axes = [{"name": nm, "length": len(ax), "origin": ax.origin, "role": ax.role,
             "induced_by": ax.induced_by} for nm, ax in ds.axes.items()]
    fields = []
    for nm, fl in ds.fields.items():
        rec = {"name": nm, "role": fl.role, "state": fl.state or None,
               "span": list(fl.span or []), "encoding": fl.encoding,
               "coverage": fl.coverage, "nullable": fl.mask is not None}
        if fl.subtype:
            rec["subtype"] = fl.subtype
        if fl.coverage == "partial":
            rec["index_axis"] = fl.index_axis
        if fl.provenance:
            rec["provenance"] = dict(fl.provenance)
        fields.append(rec)
    return {
        "source": {"path": src, "format": ff},
        "target": ({"path": dst, "format": tf} if dst else None),
        "kind": ds.kind,
        "profiles": list(ds.profiles),
        "axes": axes,
        "fields": fields,
        "dropped": list(ds.dropped),
    }


def format_report_text(rep: dict) -> str:
    """Render :func:`build_report` output as an aligned human report."""
    out = ["─ lstar conversion report ─"]
    s, t = rep["source"], rep["target"]
    line = f"source: {os.path.basename(s['path'])} ({s['format']})"
    if t:
        line += f"   →   target: {os.path.basename(t['path'])} ({t['format']})"
    out += [line, f"kind: {rep['kind']}   profiles: {', '.join(rep['profiles']) or '-'}", ""]

    aw = max((len(a["name"]) for a in rep["axes"]), default=4)
    out.append(f"axes ({len(rep['axes'])}):")
    for a in rep["axes"]:
        tag = f"   ← induced_by {a['induced_by']}" if a["induced_by"] else ""
        out.append(f"  {a['name']:<{aw}}  {a['length']:>8}  {a['origin']}"
                   f"{('/' + a['role']) if a['role'] else ''}{tag}")

    fw = max((len(f["name"]) for f in rep["fields"]), default=4)
    out += ["", f"fields ({len(rep['fields'])}):"]
    for f in rep["fields"]:
        span = " × ".join(f["span"])
        flags = []
        if f.get("state"):
            flags.append(f["state"])
        if f["coverage"] == "partial":
            flags.append(f"partial→{f.get('index_axis')}")
        if f["nullable"]:
            flags.append("nullable")
        if "provenance" in f:
            flags.append("prov{" + ",".join(list(f["provenance"])[:3]) + "}")
        tail = ("   " + " ".join(flags)) if flags else ""
        out.append(f"  {f['name']:<{fw}}  {(f['role'] or '?'):<9}  [{span}]  {f['encoding']}{tail}")

    dl = rep["dropped"]
    out += [""]
    if dl:
        out.append(f"dropped (not representable in L*): {len(dl)}")
        out += [f"  - {x}" for x in dl]
    else:
        out.append("dropped (not representable in L*): none")

    nc = rep.get("native_check")
    if nc:
        out += ["", f"native-acceptance ({nc['format']}): {nc['status']}", f"  {nc['detail']}"]
    return "\n".join(out)


def main(argv=None) -> int:
    p = argparse.ArgumentParser(
        prog="lstar", description="Convert between single-cell formats via the L* store.")
    sub = p.add_subparsers(dest="cmd", required=True)

    c = sub.add_parser("convert", help="convert SRC to DST (format detected from each path)")
    c.add_argument("src")
    c.add_argument("dst")
    c.add_argument("--from", dest="from_fmt", default=None, help="override the source format")
    c.add_argument("--to", dest="to_fmt", default=None, help="override the target format")
    c.add_argument("--backend", choices=("auto", "native", "direct"), default="auto",
                   help="native (domain package), direct (lstar's package-free codec), or auto [default]")
    c.add_argument("--report", action="store_true", help="print the full fidelity report")
    c.add_argument("--report-json", dest="report_json", metavar="FILE", default=None,
                   help="write the fidelity report as JSON to FILE")
    c.add_argument("--check", dest="check", action="store_true", default=True,
                   help="native-acceptance check on the target (open + canonical-ops smoke) [default]")
    c.add_argument("--no-check", dest="check", action="store_false",
                   help="skip the native-acceptance check")
    c.add_argument("--strict", action="store_true",
                   help="exit non-zero if the native-acceptance check fails")
    c.add_argument("-q", "--quiet", action="store_true", help="suppress the summary")

    i = sub.add_parser("inspect", help="read SRC and report its L* structure (no write)")
    i.add_argument("src")
    i.add_argument("--from", dest="from_fmt", default=None, help="override the source format")
    i.add_argument("--backend", choices=("auto", "native", "direct"), default="auto",
                   help="native, direct (package-free), or auto [default]")
    i.add_argument("--report-json", dest="report_json", metavar="FILE", default=None,
                   help="write the report as JSON to FILE")

    args = p.parse_args(argv)
    try:
        if args.cmd == "convert":
            ds, ff, tf = convert(args.src, args.dst, args.from_fmt, args.to_fmt, args.backend)
            rep = build_report(ds, args.src, ff, args.dst, tf)
            nc = None
            if args.check:
                from ._native_check import check as _native_check
                nc = _native_check(args.dst, tf)
                rep["native_check"] = nc
            if args.report_json:
                _dump_json(rep, args.report_json)
            if args.report:
                print(format_report_text(rep))
            elif not args.quiet:
                _print_summary(ds, args.src, ff, args.dst, tf)
            if nc and not (args.quiet or args.report):
                print(f"  native-acceptance ({nc['format']}): {nc['status']} — {nc['detail']}")
            if args.strict and nc and nc["status"] == "fail":
                print("lstar convert: native-acceptance check FAILED (--strict)", file=sys.stderr)
                return 3
        elif args.cmd == "inspect":
            ff = detect_format(args.src, args.from_fmt)
            ds = _read_dataset(args.src, ff, args.backend)
            rep = build_report(ds, args.src, ff, None, None)
            if args.report_json:
                _dump_json(rep, args.report_json)
            print(format_report_text(rep))
    except ConvertError as e:
        print(f"lstar {args.cmd}: error: {e}", file=sys.stderr)
        return 2
    return 0


def _dump_json(rep: dict, path: str) -> None:
    with open(path, "w") as fh:
        json.dump(rep, fh, indent=2, default=str)
