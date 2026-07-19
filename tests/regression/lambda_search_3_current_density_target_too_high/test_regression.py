"""Regression test for `examples/lambda_search_3_current_density_target_too_high`
(`target_option='max_K'`, `target_value=2e8`).

Same plasma/coil/potential setup as `lambda_search_1`; here `target_value` is
above the maximum achievable `max_K` (the lambda=0 limit, ~1.26e8 -- see
`lambda_search_1`'s test), so the legacy solver exits with `exit_code=-3`
("target too high"). `Regcoil.solve_for_target` raises `ValueError` for the
same condition.
"""

import pytest

from regcoil import CoilSurface, PlasmaSurface, Regcoil

from ..testsCommon import EQUILIBRIA


def test_target_too_high_raises():
    plasma = PlasmaSurface.from_vmec(str(EQUILIBRIA / "wout_d23p4_tm.nc"), ntheta=64, nzeta=64, mesh="full")
    coil = CoilSurface.from_uniform_offset(plasma, separation=0.5, ntheta=64, nzeta=64)
    prob = Regcoil(
        plasma, coil, mpol_potential=12, ntor_potential=12,
        net_poloidal_current=plasma.net_poloidal_current_Amperes,
        net_toroidal_current=0.0, symmetry="stellarator_symmetric",
    )

    with pytest.raises(ValueError, match="not achievable"):
        prob.solve_for_target("max_K", 2.0e8)
