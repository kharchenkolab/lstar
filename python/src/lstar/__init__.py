"""lstar — a uniform model and Zarr interchange format for single-cell/spatial omics.

A dataset is a set of Axes (entities you index by) and Fields (typed data over axes).
See ../../misc/Lstar_proposal.md for the model.
"""
from .model import Dataset, Axis, Field, Categorical, OBSERVED, DERIVED
from .zarr_io import read, write
from .validate import validate
from .lazy import LazyDense, LazyCSX, stream_col_stats
from ._engine import has_accel, show_config
from .kernels import col_sum_by_group
from .de import markers, de_bundle, de_factors, pseudobulk
from .profiles.anndata import (read_anndata, write_anndata, convert_anndata,
                               write_anndata_streamed, convert_to_h5ad)

__all__ = ["Dataset", "Axis", "Field", "Categorical", "read", "write", "validate", "OBSERVED", "DERIVED",
           "LazyDense", "LazyCSX", "stream_col_stats", "has_accel", "show_config",
           "col_sum_by_group", "markers", "de_bundle", "de_factors", "pseudobulk",
           "read_anndata", "write_anndata", "convert_anndata",
           "write_anndata_streamed", "convert_to_h5ad"]
__version__ = "0.0.1"
