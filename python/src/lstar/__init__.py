"""lstar — a uniform model and Zarr interchange format for single-cell/spatial omics.

A dataset is a set of Axes (entities you index by) and Fields (typed data over axes).
See ../../misc/Lstar_proposal.md for the model.
"""
from .model import Dataset, Axis, Field, OBSERVED, DERIVED
from .zarr_io import read, write
from .validate import validate
from .lazy import LazyDense, LazyCSX, stream_col_stats
from ._engine import has_accel, show_config
from .profiles.anndata import read_anndata, write_anndata
from .profiles.viewer import write_viewer

__all__ = ["Dataset", "Axis", "Field", "read", "write", "validate", "OBSERVED", "DERIVED",
           "LazyDense", "LazyCSX", "stream_col_stats", "has_accel", "show_config",
           "read_anndata", "write_anndata", "write_viewer"]
__version__ = "0.0.1"
