"""Regression test for `examples/regcoilPaper_figure3_NF4`
(`general_option=1`, `nlambda=1` -> a single `lambda=0` solve; full
resolution -- `ntheta_plasma=128`, so this is `@pytest.mark.slow`).

Golden values in `_golden.py` are extracted programmatically from the
example's `tests.py`.
"""

import numpy as np
import pytest

from regcoil import CoilSurface, PlasmaSurface, Regcoil

from ..tests_common import EQUILIBRIA
from ._golden import F_B, F_K, MAX_BNORMAL, MAX_K, SINGLE_VALUED_CURRENT_POTENTIAL_MN


@pytest.mark.slow
def test_li383_single_lambda():
    plasma = PlasmaSurface.from_wout(str(EQUILIBRIA / "wout_li383_1.4m.nc"), ntheta=128, nzeta=128)
    plasma.set_bnormal_from_bnorm_file(str(EQUILIBRIA / "bnorm.li383_1.4m"))
    coil = CoilSurface.from_nescin(
        str(EQUILIBRIA / "nescin.li383_realWindingSurface"), nfp=plasma.nfp, ntheta=128, nzeta=128,
    )

    prob = Regcoil(
        plasma, coil, mpol_potential=4, ntor_potential=4,
    )

    sol = prob.solve(0.0)
    np.testing.assert_allclose(sol.f_B, F_B[0], rtol=1e-4)
    np.testing.assert_allclose(sol.f_K, F_K[0], rtol=1e-4)
    np.testing.assert_allclose(sol.max_Bnormal, MAX_BNORMAL[0], rtol=1e-4)
    np.testing.assert_allclose(sol.max_K, MAX_K[0], rtol=1e-4)
    np.testing.assert_allclose(sol.solution[: len(SINGLE_VALUED_CURRENT_POTENTIAL_MN)], SINGLE_VALUED_CURRENT_POTENTIAL_MN, rtol=1e-4)
