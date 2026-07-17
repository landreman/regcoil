"""Smoke tests for the regcoil package (Phase 3 scaffold)."""

import importlib


def test_import():
    """Package can be imported."""
    regcoil = importlib.import_module("regcoil")
    assert regcoil.__version__ == "0.1.0"
