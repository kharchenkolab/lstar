"""lstar.view is a soft delegate to the pagoda3 viewer (one-way optional dependency)."""
import builtins

import pytest

import lstar


def test_view_is_exported():
    assert callable(lstar.view)
    assert "view" in lstar.__all__


def test_view_errors_with_install_hint_when_pagoda3_absent(monkeypatch):
    # Force the "pagoda3 not installed" path regardless of the environment by making its import fail.
    real_import = builtins.__import__

    def fake_import(name, *args, **kwargs):
        if name == "pagoda3" or name.startswith("pagoda3."):
            raise ImportError("No module named 'pagoda3'")
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", fake_import)
    with pytest.raises(ImportError, match="pip install pagoda3"):
        lstar.view("x.lstar.zarr")


def test_view_forwards_to_pagoda3_when_present(monkeypatch):
    # Stand in a fake pagoda3 module; view() must forward obj + kwargs and return its result.
    import sys
    import types

    seen = {}

    def fake_view(obj, **kwargs):
        seen["obj"] = obj
        seen["kwargs"] = kwargs
        return "http://viewer/url"

    fake_pagoda3 = types.ModuleType("pagoda3")
    fake_pagoda3.view = fake_view
    monkeypatch.setitem(sys.modules, "pagoda3", fake_pagoda3)

    out = lstar.view("s.lstar.zarr", prepare=False, port=1234)
    assert out == "http://viewer/url"
    assert seen["obj"] == "s.lstar.zarr"
    assert seen["kwargs"] == {"prepare": False, "port": 1234}
