"""Regression test for `examples/lambda_search_1` (`general_option=5`,
`target_option='max_K'`, `target_value=8e6`).

The legacy Brent search visits lambda in a different order than
`Regcoil.solve_for_target` (a direct bisection on the closed-form family,
ADR-021, see PHASES.md Phase 8), so only the converged endpoint and the
lambda=0 / lambda=inf limits are checked here, not the intermediate
sequence. Golden values are read by hand from `examples/lambda_search_1/tests.py`.
"""

import numpy as np

from regcoil import CoilSurface, PlasmaSurface, Regcoil

from ..tests_common import EQUILIBRIA


def _build_problem():
    plasma = PlasmaSurface.from_wout(str(EQUILIBRIA / "wout_d23p4_tm.nc"), ntheta=64, nzeta=64)
    coil = CoilSurface.from_uniform_offset(
        plasma, separation=0.5, ntheta=64, nzeta=64, standard_toroidal_angle=True,
        theta_reparameterization=None,
    )
    return Regcoil(
        plasma, coil, mpol_potential=12, ntor_potential=12,
    )


def test_lambda_zero_and_infinite_limits():
    # lambda=1e200 (lambda=inf) and lambda=0 -- the first two entries of the
    # legacy lambda array for this example.
    prob = _build_problem()

    sol_inf = prob.solve(np.inf)
    np.testing.assert_allclose(sol_inf.f_B, 13.1506816246147, rtol=0.001)
    np.testing.assert_allclose(sol_inf.f_K, 1.21300287142288e15, rtol=0.001)
    np.testing.assert_allclose(sol_inf.max_Bnormal, 0.800268516799401, rtol=0.001)
    np.testing.assert_allclose(sol_inf.max_K, 2754022.92888632, rtol=0.001)

    sol_0 = prob.solve(0.0)
    np.testing.assert_allclose(sol_0.f_B, 7.71670836033572e-06, rtol=0.001)
    np.testing.assert_allclose(sol_0.f_K, 1.99560519855456e16, rtol=0.001)
    np.testing.assert_allclose(sol_0.max_Bnormal, 0.00131646260910012, rtol=0.001)
    np.testing.assert_allclose(sol_0.max_K, 126095244.116546, rtol=0.001)


def test_solve_for_target_max_K():
    prob = _build_problem()
    sol = prob.solve_for_target("max_K", 8.0e6)
    np.testing.assert_allclose(sol.max_K, 8.0e6, rtol=1e-3)
    np.testing.assert_allclose(sol.f_B, 0.337712583344148, rtol=1e-4)
