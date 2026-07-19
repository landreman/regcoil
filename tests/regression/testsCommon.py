"""Shared helpers for tests/regression/*/test_regression.py (Phase 8).

These tests build the problem directly via the Python object model
(`regcoil.PlasmaSurface`/`CoilSurface`/`Regcoil`) and compare against golden
values taken by hand from the legacy Fortran solver -- the same reference
values recorded in the corresponding `examples/*/tests.py` (see PHASES.md
Phase 8). There is no more legacy-executable subprocess or netCDF output to
read.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

EQUILIBRIA = Path(__file__).resolve().parents[2] / "equilibria"


def legacy_lambda_array(nlambda, lambda_min, lambda_max):
    """Port of the legacy fixed lambda scan (`regcoil_compute_lambda.f90`,
    `general_option=1`): `lambda[0] = 0`, then `nlambda-1` values log-spaced
    from `lambda_min` to `lambda_max`.
    """
    if nlambda == 1:
        return np.array([0.0])
    lambdas = np.empty(nlambda)
    lambdas[0] = 0.0
    j = np.arange(1, nlambda)
    lambdas[1:] = lambda_min * np.exp(np.log(lambda_max / lambda_min) * (j - 1) / (nlambda - 2))
    return lambdas
