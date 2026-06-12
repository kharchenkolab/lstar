#!/usr/bin/env python3
"""Streaming vs in-memory conversion: peak memory and wall time as the dataset grows.

Subsamples a real .h5ad to several sizes and, for each, converts it to an L* store two ways --
**eager** (load the whole AnnData, then write) and **streaming** (`convert_anndata`, backed read +
block-by-block write) -- measuring peak RSS (via `/usr/bin/time -v`, in a fresh subprocess so the
off-heap C++/zarr allocations are counted) and wall time. Produces a 2-panel figure: peak memory and
wall time vs dataset size. The point: streaming holds memory ~flat as the data grows, for a modest
time cost.

Usage:
  python3 examples/streaming_scaling.py [source.h5ad] [out.png]
  python3 examples/streaming_scaling.py --worker {eager|stream} <h5ad> <store>   # internal
"""
import os
import sys
import time
import shutil
import subprocess

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "python", "src"))

DEFAULT_SRC = "/tmp/tms_droplet_full.h5ad"   # Tabula Muris Senis droplet atlas (~245k cells)
SIZE_FRACS = [0.1, 0.25, 0.5, 1.0]           # subsample to these fractions of cells; nnz scales ~linearly
CHUNK_ELEMS = 2_000_000


def worker(mode, h5, store):
    """One conversion in this (isolated) process; prints WALL=<seconds>."""
    import anndata as ad
    import lstar
    from lstar.profiles.anndata import read_anndata, convert_anndata
    if os.path.exists(store):
        shutil.rmtree(store)
    t = time.time()
    if mode == "eager":
        lstar.write(read_anndata(ad.read_h5ad(h5)), store)
    else:
        convert_anndata(h5, store, chunk_elems=CHUNK_ELEMS)
    print("WALL=%.2f" % (time.time() - t))


def measure(mode, h5, store):
    """Run the worker under /usr/bin/time -v; return (wall_s, peak_rss_mb)."""
    p = subprocess.run(["/usr/bin/time", "-v", sys.executable, os.path.abspath(__file__),
                        "--worker", mode, h5, store],
                       capture_output=True, text=True)
    out = p.stdout + p.stderr
    wall = next((float(l.split("=")[1]) for l in out.splitlines() if l.startswith("WALL=")), None)
    rss = next((int(l.split(":")[1]) / 1024.0 for l in out.splitlines() if "Maximum resident" in l), None)
    if wall is None or rss is None:
        raise RuntimeError("measurement failed for %s:\n%s" % (mode, out[-800:]))
    return wall, rss


def make_subsamples(src):
    """Read the source once (backed) and write a minimal h5ad (X only) at each size fraction."""
    import anndata as ad
    import scipy.sparse as sp
    a = ad.read_h5ad(src, backed="r")
    ntot = a.n_obs
    out = []  # (cells, nnz, h5_path)
    for frac in SIZE_FRACS:
        n = int(round(ntot * frac))
        p = "/tmp/scale_%d.h5ad" % n
        if os.path.exists(p):                    # reuse a prepared subsample (re-runs are fast)
            nnz = int(sp.csr_matrix(ad.read_h5ad(p).X).nnz)
        else:
            Xn = sp.csr_matrix(a.X[:n])          # backed slice -> reads only n rows from disk
            Xn.data = Xn.data.astype("f4")
            ad.AnnData(X=Xn).write_h5ad(p)
            nnz = int(Xn.nnz)
        out.append((n, nnz, p))
        print("  prepared %d cells (%d nnz) -> %s" % (n, nnz, p))
    try:
        a.file.close()
    except Exception:
        pass
    return out


def main():
    src = sys.argv[1] if len(sys.argv) > 1 and not sys.argv[1].startswith("-") else DEFAULT_SRC
    out_png = sys.argv[2] if len(sys.argv) > 2 else os.path.join(HERE, "..", "docs", "img", "streaming_scaling.png")
    rows = []  # (cells, nnz, mode, wall, rss)
    for n, nnz, h5 in make_subsamples(src):
        for mode in ("eager", "stream"):
            wall, rss = measure(mode, h5, "/tmp/scale_%s.lstar.zarr" % mode)
            rows.append((n, nnz, mode, wall, rss))
            print("cells=%-7d nnz=%-11d %-7s wall=%6.1fs peakRSS=%7.0f MB" % (n, nnz, mode, wall, rss))
    plot(rows, out_png)
    print("\nwrote %s" % os.path.abspath(out_png))


def plot(rows, out_png):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    xs = sorted({r[0] for r in rows})
    xk = [x / 1000.0 for x in xs]
    tot_nnz = max(r[1] for r in rows)
    def series(mode, idx):
        return [next(r[idx] for r in rows if r[0] == x and r[2] == mode) for x in xs]
    style = {"eager": dict(color="#c0392b", marker="o", label="in-memory (eager)"),
             "stream": dict(color="#2471a3", marker="s", label="streaming")}
    fig, (axm, axt) = plt.subplots(1, 2, figsize=(8.4, 3.4))
    for mode in ("eager", "stream"):
        axm.plot(xk, series(mode, 4), lw=2, ms=6, **style[mode])
        axt.plot(xk, series(mode, 3), lw=2, ms=6, **style[mode])
    axm.set(xlabel="cells (thousands)", ylabel="peak memory (MB)", title="Peak memory")
    axt.set(xlabel="cells (thousands)", ylabel="wall time (s)", title="Wall time")
    for ax in (axm, axt):
        ax.grid(alpha=0.3); ax.set_xlim(0, max(xk) * 1.05); ax.set_ylim(bottom=0)
    axm.legend(frameon=False, fontsize=9, loc="upper left")
    fig.suptitle("Streaming vs in-memory: h5ad → L* conversion (Tabula Muris Senis droplet, up to %.0fM nonzeros)"
                 % (tot_nnz / 1e6), fontsize=10)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    os.makedirs(os.path.dirname(out_png), exist_ok=True)
    fig.savefig(out_png, dpi=130)


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "--worker":
        worker(sys.argv[2], sys.argv[3], sys.argv[4])
    else:
        main()
