"""Smoke test: package is importable after install / editable layout."""

import regcoil


def test_version():
    assert regcoil.__version__
