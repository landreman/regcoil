"""Regression test for
`examples/regcoilPaper_figure10d_but_geometry_option_coil_4_loRes`
(`geometry_option_coil=4`, uniform offset with the legacy constant-arclength
theta reparametrization).

Since ADR-031 this is the real `geometry_option_coil=4` construction:
`standard_toroidal_angle=True` (the root-solved toroidal angle, legacy `=2`)
plus `theta_reparameterization="uniform_arclength"` (the constant-arclength
theta the legacy code reached by fixed-point iteration, reached here by a
single quadrature).

The golden values are the legacy example's own, which are shared with
`geometry_option_coil=2` and the constant-arclength-nescin variant: they are
numerically identical because *that is the example's point* -- the physics
does not depend on the coil-surface parametrization. So this test is the
acceptance test for the reparameterization: it checks both that the legacy
`=4` numbers are reproduced, and (via the shared goldens) that reparameterizing
did not move the physical surface.
"""

import numpy as np

from regcoil import CoilSurface, PlasmaSurface, Regcoil

from ..tests_common import EQUILIBRIA, lambda_array


def test_uniform_offset_lores_matches_constant_arclength_angle():
    plasma = PlasmaSurface.from_wout(str(EQUILIBRIA / "wout_d23p4_tm.nc"), ntheta=64, nzeta=64)
    plasma.set_bnormal_from_bnorm_file(str(EQUILIBRIA / "bnorm.d23p4_tm"))
    coil = CoilSurface.from_uniform_offset(
        plasma, separation=0.5, ntheta=64, nzeta=64, standard_toroidal_angle=True,
        theta_reparameterization="uniform_arclength",
    )
    # The legacy `constant_arclength_tolerance` thresholded exactly this.
    assert coil.theta_map.diagnostics["max_dl_spread"] < 1e-4

    prob = Regcoil(
        plasma, coil, mpol_potential=12, ntor_potential=12,
    )

    lambdas = lambda_array(nlambda=10, lambda_min=1e-15, lambda_max=1e-14)
    sols = prob.scan(lambdas)

    f_B = np.array([sol.f_B for sol in sols])
    np.testing.assert_allclose(
        f_B[1:],
        [0.174519878306313, 0.233148551834137, 0.310446213587783, 0.411565852326815,
         0.542616084841207, 0.710568315783369, 0.922890554110523, 1.18673864681126,
         1.50770781013225],
        rtol=0.03,
    )

    f_K = np.array([sol.f_K for sol in sols])
    np.testing.assert_allclose(
        f_K[1:],
        [1.74957088182873e15, 1.6989658516062e15, 1.64892520707378e15, 1.59982527129183e15,
         1.55209554042516e15, 1.50621123410416e15, 1.46269693947801e15, 1.42212837760326e15,
         1.38509938973785e15],
        rtol=0.01,
    )

    max_Bnormal = np.array([sol.max_Bnormal for sol in sols])
    np.testing.assert_allclose(
        max_Bnormal[1:],
        [0.117555544291034, 0.134561224876993, 0.153923792979921, 0.176347144447289,
         0.201603899107085, 0.229977435728317, 0.261839981656698, 0.296625161466877,
         0.333988791941713],
        rtol=0.03,
    )

    max_K = np.array([sol.max_K for sol in sols])
    np.testing.assert_allclose(
        max_K[1:],
        [8824481.60102707, 8355450.23653953, 7884212.39105115, 7412969.17310808,
         6945281.57578389, 6485916.75531679, 6040435.92861007, 5614638.75966927,
         5214018.48872153],
        rtol=0.06,
    )
