"""Regression test for `examples/axisymmetrySanityTest_chi2K_regularization`.

Both plasma and coil surfaces are axisymmetric circular tori (different major
radius), so the single-valued part of the current potential and B_normal
should vanish to near machine precision at every lambda -- a strong,
resolution-independent sanity check of the basis-function/matrix assembly
that does not depend on any golden legacy value.
"""

import numpy as np

from regcoil import CoilSurface, PlasmaSurface, Regcoil

from ..tests_common import lambda_array


def test_axisymmetric_solution_vanishes():
    plasma = PlasmaSurface.circular_torus(R0=6, a=2, nfp=3, ntheta=32, nzeta=32)
    plasma.net_poloidal_current = 1.0e6
    coil = CoilSurface.circular_torus(R0=6.5, a=4, nfp=3, ntheta=32, nzeta=32)

    prob = Regcoil(
        plasma, coil, mpol_potential=6, ntor_potential=7,
    )
    assert prob.nbf == 97  # legacy: single_valued_current_potential_mn has length 97

    lambdas = lambda_array(nlambda=3, lambda_min=1e-15, lambda_max=1e100)
    for sol in prob.scan(lambdas):
        assert abs(sol.chi2_B) < 1e-10
        assert abs(sol.max_Bnormal) < 1e-10
        np.testing.assert_allclose(sol.solution, 0, atol=3e-10)
