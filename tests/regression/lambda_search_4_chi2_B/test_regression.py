"""Regression test for `examples/lambda_search_4_chi2_B`
(`general_option=5`, `target_option='chi2_B'`, `target_value=0.1`).

See `lambda_search_1` for why only the converged endpoint is checked, not
the legacy Brent search's intermediate lambda sequence.
"""

import numpy as np

from regcoil import CoilSurface, PlasmaSurface, Regcoil

from ..testsCommon import EQUILIBRIA


def test_solve_for_target_chi2_B():
    plasma = PlasmaSurface.from_vmec(str(EQUILIBRIA / "wout_d23p4_tm.nc"), ntheta=64, nzeta=64, mesh="full")
    coil = CoilSurface.from_uniform_offset(plasma, separation=0.5, ntheta=64, nzeta=64, standard_toroidal_angle=True)
    prob = Regcoil(
        plasma, coil, mpol_potential=12, ntor_potential=12,
    )

    sol = prob.solve_for_target("chi2_B", 0.1)
    np.testing.assert_allclose(sol.chi2_B, 0.0999995031486614, rtol=1e-4)
