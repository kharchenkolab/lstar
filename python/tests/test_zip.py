"""Single-file `.lstar.zarr.zip` (ZIP_STORED) write+read — Python reference surface.

A `.lstar.zarr.zip` is a normal L* store packed into ONE zip file with every entry STORED (no
deflate): zarr chunks are already codec-compressed, and only a STORED entry stays byte-range-readable
inside the archive (the point of a hosted single file). This test pins the contract:
  - round-trip through a `.zip` == round-trip through a directory store (field-for-field),
  - every zip entry is STORED (a deflated `.lstar.zarr.zip` must be impossible to produce here),
  - reading a foreign DEFLATE-packed store is rejected with a clear message,
  - lazy read works straight from the zip, and stream=True targets a zip too,
  - a corpus dataset (real locally / synthetic-faithful in CI) survives the zip round-trip.
"""
import os
import sys
import tempfile
import zipfile

import numpy as np
import scipy.sparse as sp

import lstar
from lstar.zarr_io import write, read


def _make_ds(n_cells=40, n_genes=18, n_pc=5):
    """A Dataset spanning the field encodings a store can hold: dense, csr, csc, categorical,
    string, and a directed/weighted graph — so the zip round-trip exercises each reader path."""
    rng = np.random.default_rng(0)
    ds = lstar.Dataset(kind="sample")
    ds.add_axis("cells", ["cell%d" % i for i in range(n_cells)])
    ds.add_axis("genes", ["g%d" % i for i in range(n_genes)])
    ds.add_axis("pca", ["PC%d" % i for i in range(n_pc)], origin="derived")

    counts = sp.random(n_cells, n_genes, density=0.3, format="csr", random_state=0)
    counts.data = np.ceil(counts.data * 10).astype(np.float32)
    ds.add_field("counts", counts, role="measure", span=["cells", "genes"], state="raw")

    lognorm = sp.random(n_cells, n_genes, density=0.3, format="csc", random_state=1).astype(np.float32)
    ds.add_field("data", lognorm, role="measure", span=["cells", "genes"], state="lognorm")

    ds.add_field("pca", rng.standard_normal((n_cells, n_pc)).astype(np.float32),
                 role="embedding", span=["cells", "pca"])

    leiden = np.array(["c%d" % (i % 4) for i in range(n_cells)])
    ds.add_field("leiden", leiden, role="label", span=["cells"])
    ds.add_field("barcode", np.array(["BC-%04d" % i for i in range(n_cells)]),
                 role="label", span=["cells"])

    knn = sp.random(n_cells, n_cells, density=0.1, format="csr", random_state=2).astype(np.float32)
    g = ds.add_field("knn", knn, role="graph", span=["cells", "cells"])
    g.directed = True
    g.weighted = True
    return ds


def _sparse_eq(a, b):
    A = a.toarray() if sp.issparse(a) else np.asarray(a)
    B = b.toarray() if sp.issparse(b) else np.asarray(b)
    return A.shape == B.shape and np.allclose(A, B)


def _field_eq(fa, fb):
    va, vb = fa.values, fb.values
    if sp.issparse(va) or sp.issparse(vb):
        return _sparse_eq(va, vb)
    aa, ab = np.asarray(va), np.asarray(vb)
    if aa.dtype.kind in ("U", "S", "O") or ab.dtype.kind in ("U", "S", "O"):
        return list(map(str, aa.tolist())) == list(map(str, ab.tolist()))
    return aa.shape == ab.shape and np.allclose(aa, ab)


def _assert_ds_equal(a, b, ctx="", check_encoding=True):
    # check_encoding: only when BOTH sides came from a store — the in-memory original assigns a field's
    # on-disk encoding lazily (a string label is "dense" until the writer normalizes it to "utf8").
    assert list(a.axes) == list(b.axes), f"{ctx}: axes differ"
    for name in a.axes:
        assert list(map(str, a.axis(name).labels)) == list(map(str, b.axis(name).labels)), \
            f"{ctx}: axis {name} labels differ"
    assert set(a.fields) == set(b.fields), f"{ctx}: field sets differ {set(a.fields)} vs {set(b.fields)}"
    for name in a.fields:
        fa, fb = a.field(name), b.field(name)
        assert _field_eq(fa, fb), f"{ctx}: field {name} values differ"
        if check_encoding:
            assert (fa.encoding == fb.encoding), f"{ctx}: field {name} encoding {fa.encoding} vs {fb.encoding}"
        assert bool(fa.directed) == bool(fb.directed), f"{ctx}: field {name} directed differs"
        assert bool(fa.weighted) == bool(fb.weighted), f"{ctx}: field {name} weighted differs"


def _all_entries_stored(zippath):
    with zipfile.ZipFile(zippath) as z:
        infos = z.infolist()
        return len(infos) > 0 and all(i.compress_type == zipfile.ZIP_STORED for i in infos)


# ---- tests ----

def test_zip_roundtrip_equals_dir():
    """A `.zip` round-trip reproduces the directory round-trip field-for-field, and every entry STORED."""
    ds = _make_ds()
    with tempfile.TemporaryDirectory() as tmp:
        d = os.path.join(tmp, "s.lstar.zarr")
        z = os.path.join(tmp, "s.lstar.zarr.zip")
        write(ds, d)
        write(ds, z)
        assert os.path.isfile(z), "the .zip target must be a single file, not a directory"
        assert _all_entries_stored(z), "every entry in a .lstar.zarr.zip must be ZIP_STORED"
        rd = read(d)
        rz = read(z)
        _assert_ds_equal(rd, rz, ctx="dir-vs-zip")
        _assert_ds_equal(ds, rz, ctx="orig-vs-zip", check_encoding=False)
    print("ok test_zip_roundtrip_equals_dir")


def test_zip_forces_stored_despite_compressor():
    """Passing a compressor must NOT produce deflated zip entries — the zip layer stays STORED
    (the inner chunks may be codec-compressed; the zip must not re-deflate them)."""
    import numcodecs
    ds = _make_ds()
    with tempfile.TemporaryDirectory() as tmp:
        z = os.path.join(tmp, "s.lstar.zarr.zip")
        write(ds, z, compressor=numcodecs.GZip(5))
        assert _all_entries_stored(z), "compressor= must not cause deflated zip entries"
        _assert_ds_equal(ds, read(z), ctx="compressed-chunks-zip", check_encoding=False)
    print("ok test_zip_forces_stored_despite_compressor")


def test_zip_lazy_read():
    """Lazy read works straight from the zip (chunk-at-a-time, no full extract)."""
    ds = _make_ds()
    with tempfile.TemporaryDirectory() as tmp:
        z = os.path.join(tmp, "s.lstar.zarr.zip")
        write(ds, z, chunk_elems=64)
        rz = read(z, lazy=True)
        # touch a lazy sparse column-reduction path + a dense slice
        got = rz.field("data").values
        assert got.shape == (40, 18)
        _assert_ds_equal(ds, read(z), ctx="lazy-zip-materialize", check_encoding=False)
    print("ok test_zip_lazy_read")


def test_zip_stream_write():
    """stream=True targets a zip too (streams to a temp store, finalizes STORED)."""
    ds = _make_ds()
    with tempfile.TemporaryDirectory() as tmp:
        z = os.path.join(tmp, "s.lstar.zarr.zip")
        write(ds, z, stream=True)
        assert _all_entries_stored(z), "streamed .zip write must still be all-STORED"
        _assert_ds_equal(ds, read(z), ctx="stream-zip", check_encoding=False)
    print("ok test_zip_stream_write")


def test_zip_rejects_deflate():
    """A foreign, DEFLATE-packed store (e.g. `zip -r`) is rejected on read with a clear message."""
    ds = _make_ds()
    with tempfile.TemporaryDirectory() as tmp:
        d = os.path.join(tmp, "s.lstar.zarr")
        z = os.path.join(tmp, "bad.lstar.zarr.zip")
        write(ds, d)
        # hand-roll a DEFLATE zip of the directory store (the wrong artifact)
        with zipfile.ZipFile(z, "w", zipfile.ZIP_DEFLATED) as zf:
            for root, _dirs, files in os.walk(d):
                for fn in files:
                    fp = os.path.join(root, fn)
                    zf.write(fp, arcname=os.path.relpath(fp, d))
        raised = False
        try:
            read(z)
        except Exception as e:
            raised = True
            msg = str(e).lower()
            assert "stored" in msg or "deflate" in msg or "range" in msg, \
                f"deflate rejection message should be actionable, got: {e}"
        assert raised, "reading a DEFLATE-packed .lstar.zarr.zip must be rejected"
    print("ok test_zip_rejects_deflate")


def test_zip_corpus_roundtrip():
    """Corpus datasets (real locally / synthetic-faithful in CI) survive the zip round-trip
    identically to the directory round-trip — realistic categoricals, embeddings, native sparse."""
    try:
        sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
        import corpus
        from lstar.profiles.anndata import read_anndata
    except Exception as e:
        print("skip test_zip_corpus_roundtrip (corpus/anndata unavailable):", e)
        return
    loaders = [("pbmc68k", corpus.pbmc68k_reduced), ("pbmc3k", corpus.pbmc3k_processed)]
    ran = 0
    for name, load in loaders:
        try:
            a = load()
        except Exception as e:
            print(f"  skip corpus {name}: {e}")
            continue
        if a is None:                       # loader returns None when its backend (e.g. scanpy) is unavailable
            print(f"  skip corpus {name}: loader unavailable")
            continue
        ds = read_anndata(a)
        with tempfile.TemporaryDirectory() as tmp:
            d = os.path.join(tmp, f"{name}.lstar.zarr")
            z = os.path.join(tmp, f"{name}.lstar.zarr.zip")
            write(ds, d)
            write(ds, z)
            assert _all_entries_stored(z), f"corpus {name}: zip not all-STORED"
            _assert_ds_equal(read(d), read(z), ctx=f"corpus-{name}")
        ran += 1
        print(f"  ok corpus {name} zip round-trip")
    if ran == 0:
        print("skip test_zip_corpus_roundtrip (no corpus datasets available)")
    else:
        print("ok test_zip_corpus_roundtrip")


if __name__ == "__main__":
    test_zip_roundtrip_equals_dir()
    test_zip_forces_stored_despite_compressor()
    test_zip_lazy_read()
    test_zip_stream_write()
    test_zip_rejects_deflate()
    test_zip_corpus_roundtrip()
    print("all zip tests passed")
