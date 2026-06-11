"""Version recognition: profiles detect the source format version and degrade gracefully.

External formats have many versions with different object layouts. The profiles must recognize
which one they're reading (recorded in ds.profiles) and adapt to layout variants -- here, the
AnnData `.raw` slot, which older pipelines use to stash pre-HVG raw counts over a *different*
gene set. It must get its own axis and round-trip, not be forced onto `genes` or dropped.
"""
import os
import tempfile

import numpy as np
import pandas as pd
import scipy.sparse as sp

from lstar import read, write, validate
from lstar.profiles.anndata import read_anndata, write_anndata


def test_records_anndata_version():
    import anndata as ad
    a = ad.AnnData(X=sp.random(10, 5, density=0.3, format="csr"))
    ds = read_anndata(a)
    vers = [p for p in ds.profiles if p.startswith("anndata@") and p != "anndata@0.1"]
    assert vers, ds.profiles                      # the detected library version is recorded
    print("records anndata version: %s" % vers[0])


def test_raw_with_divergent_genes():
    import anndata as ad
    rawX = sp.random(30, 50, density=0.2, format="csr")
    raw = ad.AnnData(X=rawX, var=pd.DataFrame(index=[f"g{i}" for i in range(50)]))
    a = ad.AnnData(X=rawX[:, :20].copy(), var=pd.DataFrame(index=[f"g{i}" for i in range(20)]))
    a.raw = raw

    ds = read_anndata(a)
    assert "genes_raw" in ds.axes and len(ds.axes["genes_raw"]) == 50
    rf = ds.fields["raw"]
    assert rf.span == ["cells", "genes_raw"] and rf.state == "raw"
    assert not validate(ds)

    p = os.path.join(tempfile.mkdtemp(), "raw.lstar.zarr")
    write(ds, p)
    a2 = write_anndata(read(p))
    assert a2.raw is not None and a2.raw.shape == (30, 50) and a2.shape == (30, 20)
    print("raw (divergent gene set) recognized + round-trips: X %s, raw %s"
          % (a2.shape, a2.raw.shape))


if __name__ == "__main__":
    test_records_anndata_version()
    test_raw_with_divergent_genes()
