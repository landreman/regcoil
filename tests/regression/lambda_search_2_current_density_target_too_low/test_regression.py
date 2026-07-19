"""Regression test for `examples/lambda_search_2_current_density_target_too_low`
(`target_option='max_K'`, `target_value=1e6`).

Same plasma/coil/potential setup as `lambda_search_1`; here `target_value` is
below the minimum achievable `max_K` (the lambda=inf limit, ~2.75e6 -- see
`lambda_search_1`'s test), so the legacy solver exits with `exit_code=-2`
("target too low"). `Regcoil.solve_for_target` raises `ValueError` for the
same condition.
"""

import pytest

from regcoil import CoilSurface, PlasmaSurface, Regcoil

from ..tests_common import EQUILIBRIA


def test_target_too_low_raises():
    plasma = PlasmaSurface.from_vmec(str(EQUILIBRIA / "wout_d23p4_tm.nc"), ntheta=64, nzeta=64, mesh="full")
    coil = CoilSurface.from_uniform_offset(plasma, separation=0.5, ntheta=64, nzeta=64, standard_toroidal_angle=True)
    prob = Regcoil(
        plasma, coil, mpol_potential=12, ntor_potential=12,
    )

    with pytest.raises(ValueError, match="not achievable"):
        prob.solve_for_target("max_K", 1.0e6)
