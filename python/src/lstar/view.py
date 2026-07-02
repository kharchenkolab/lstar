"""Soft delegate to the pagoda3 viewer.

lstar is the substrate; interactive viewing lives in the *separate* ``pagoda3`` package (which
depends on lstar, not the other way around). The dependency stays one-way: lstar never *requires*
the viewer — it lazily imports ``pagoda3`` and forwards to it when present, and otherwise raises a
clear install hint. The real logic (coerce -> prep -> serve -> open) lives in :func:`pagoda3.view`.
"""
from __future__ import annotations


def view(obj, **kwargs):
    """Open ``obj`` in the pagoda3 viewer (delegates to the ``pagoda3`` package).

    Convenience shim for an lstar session: forwards to :func:`pagoda3.view`. ``obj`` may be a
    ``*.lstar.zarr`` store path, an :class:`lstar.Dataset`, or an AnnData. Requires the separate
    ``pagoda3`` package (``pip install pagoda3``); lstar only optionally depends on it, so a plain
    lstar install has no viewer weight.

    Any keyword arguments (``prepare``, ``local``, ``host``, ``port``, ``open_browser`` …) pass
    straight through to :func:`pagoda3.view`.
    """
    try:
        from pagoda3 import view as _view
    except ImportError as e:  # keep the dependency one-way: pagoda3 is optional for lstar
        raise ImportError(
            "Interactive viewing needs the 'pagoda3' package:  pip install pagoda3"
        ) from e
    return _view(obj, **kwargs)
