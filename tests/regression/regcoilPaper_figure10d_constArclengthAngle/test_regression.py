"""Regression test for `examples/regcoilPaper_figure10d_constArclengthAngle`
(`geometry_option_coil=3`, a nescin file for the *constant-arclength*
theta reparametrization of the uniform-offset coil surface, full resolution
-- `@pytest.mark.slow`).

The point of this legacy example is that the physics (chi2_B, chi2_K,
max_Bnormal, max_K) does not depend on the coil-surface parametrization,
even though the individual Fourier mode coefficients of the current
potential do (different theta coordinate on the same physical surface) --
hence the generous chi2 tolerances and the separate golden
`single_valued_current_potential_mn` (distinct from the plain nescin file's
in the `_originalAngle` sibling).

Golden values in `_golden.py` are extracted programmatically from the
example's `tests.py` (large arrays -- not hand-transcribed).
"""

import numpy as np
import pytest

from regcoil import CoilSurface, PlasmaSurface, Regcoil

from ..testsCommon import EQUILIBRIA, legacy_lambda_array
from ._golden import CHI2_B, CHI2_K, LAMBDA, MAX_BNORMAL, MAX_K, SINGLE_VALUED_CURRENT_POTENTIAL_MN


@pytest.mark.slow
def test_nescin_const_arclength_highres():
    plasma = PlasmaSurface.from_vmec(str(EQUILIBRIA / "wout_d23p4_tm.nc"), ntheta=128, nzeta=128, mesh="full")
    plasma.set_bnormal_from_bnorm_file(str(EQUILIBRIA / "bnorm.d23p4_tm"))
    coil = CoilSurface.from_nescin(
        str(EQUILIBRIA / "nescin.d23p4_tm_uniform_0.5m_offset_eqarc"), nfp=plasma.nfp, ntheta=128, nzeta=128,
    )

    prob = Regcoil(
        plasma, coil, mpol_potential=32, ntor_potential=32,
    )

    lambdas = legacy_lambda_array(nlambda=10, lambda_min=1e-15, lambda_max=1e-14)
    np.testing.assert_allclose(lambdas, LAMBDA, rtol=1e-12)
    sols = prob.scan(lambdas)

    np.testing.assert_allclose([sol.chi2_B for sol in sols][1:], CHI2_B, rtol=3e-3)
    np.testing.assert_allclose([sol.chi2_K for sol in sols][1:], CHI2_K, rtol=1e-4)
    np.testing.assert_allclose([sol.max_Bnormal for sol in sols][1:], MAX_BNORMAL, rtol=1e-3)
    np.testing.assert_allclose([sol.max_K for sol in sols][1:], MAX_K, rtol=1e-3)
    np.testing.assert_allclose(
        sols[1].solution[: len(SINGLE_VALUED_CURRENT_POTENTIAL_MN)],
        SINGLE_VALUED_CURRENT_POTENTIAL_MN,
        rtol=1e-4,
        atol=1e-6,
    )
