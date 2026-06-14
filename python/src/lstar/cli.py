"""``lstar convert`` — one-command conversion between single-cell formats, bridged by the L* store.

``lstar convert SRC DST`` detects each format from its path, reads SRC into the L* model, and writes DST
from it. The L* dataset is the universal intermediate, so a conversion is just ``write_Y(read_X(obj))``
and what a target cannot hold is recorded in ``ds.dropped`` (visible, never silently lost) rather than
dropped on the floor.

Same-language conversions (AnnData / MuData / store, all Python-side) run in-process. Cross-language ones
(Seurat / SCE, materialized in R) are bridged by a temporary on-disk store — wired in a later step; this
build handles the Python-side formats.

Format detection (override with ``--from`` / ``--to``):
  ``.h5ad`` → anndata · ``.h5mu`` → mudata · ``.rds`` → rds (Seurat/SCE, sniffed R-side) ·
  ``.lstar.zarr`` / ``.zarr`` / a Zarr directory → store
"""
from __future__ import annotations

import argparse
import json
import os
import sys


class ConvertError(Exception):
    """A user-facing conversion error (bad path, undetectable format, unsupported route)."""


# Which language materializes each format.
_PY = {"anndata", "mudata", "store"}
_R = {"seurat", "sce", "rds"}

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


def _load_py(src: str, fmt: str):
    """Read a Python-side source into an L* :class:`Dataset`."""
    if not os.path.exists(src):
        raise ConvertError(f"source not found: {src}")
    import lstar
    if fmt == "store":
        return lstar.read(src)
    if fmt == "anndata":
        try:
            import anndata as ad
        except ImportError:
            raise ConvertError("reading AnnData (.h5ad) needs the 'anndata' package (pip install anndata)")
        from lstar.profiles.anndata import read_anndata
        return read_anndata(ad.read_h5ad(src))
    if fmt == "mudata":
        try:
            import mudata as md
        except ImportError:
            raise ConvertError("reading MuData (.h5mu) needs the 'mudata' package (pip install mudata)")
        from lstar.profiles.mudata import read_mudata
        return read_mudata(md.read_h5mu(src))
    raise ConvertError(f"{fmt!r} is not a Python-side source format")


def _emit_py(ds, dst: str, fmt: str) -> None:
    """Write an L* :class:`Dataset` out to a Python-side target."""
    import lstar
    if fmt == "store":
        lstar.write(ds, dst)
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


def convert(src: str, dst: str, from_fmt: str | None = None, to_fmt: str | None = None):
    """Convert *src* → *dst*. Returns ``(ds, from_fmt, to_fmt)`` with the bridging L* dataset."""
    ff = detect_format(src, from_fmt)
    tf = detect_format(dst, to_fmt)
    if ff in _R or tf in _R:
        raise ConvertError(
            "Seurat/SCE (R-side) conversion is not wired in this build yet — "
            "use the Python formats (anndata | mudata | store) for now")
    ds = _load_py(src, ff)
    _emit_py(ds, dst, tf)
    return ds, ff, tf


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
    c.add_argument("--report", action="store_true", help="print the full fidelity report")
    c.add_argument("--report-json", dest="report_json", metavar="FILE", default=None,
                   help="write the fidelity report as JSON to FILE")
    c.add_argument("-q", "--quiet", action="store_true", help="suppress the summary")

    i = sub.add_parser("inspect", help="read SRC and report its L* structure (no write)")
    i.add_argument("src")
    i.add_argument("--from", dest="from_fmt", default=None, help="override the source format")
    i.add_argument("--report-json", dest="report_json", metavar="FILE", default=None,
                   help="write the report as JSON to FILE")

    args = p.parse_args(argv)
    try:
        if args.cmd == "convert":
            ds, ff, tf = convert(args.src, args.dst, args.from_fmt, args.to_fmt)
            rep = build_report(ds, args.src, ff, args.dst, tf)
            if args.report_json:
                _dump_json(rep, args.report_json)
            if args.report:
                print(format_report_text(rep))
            elif not args.quiet:
                _print_summary(ds, args.src, ff, args.dst, tf)
        elif args.cmd == "inspect":
            ff = detect_format(args.src, args.from_fmt)
            if ff in _R:
                raise ConvertError("inspecting Seurat/SCE (R-side) is not wired in this build yet")
            ds = _load_py(args.src, ff)
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
