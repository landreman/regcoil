"""Unit tests for `Regcoil`/`Solution` (Phase 8): object-model behavior that
doesn't need a golden legacy comparison (see tests/regression/ for those).
"""

import numpy as np
import pytest

from regcoil import CoilSurface, PlasmaSurface, Regcoil


def _small_problem(ntheta=8, nzeta=8, mpol=3, ntor=2, R0=6.0, a_plasma=2.0, a_coil=3.0, nfp=3):
    plasma = PlasmaSurface.circular_torus(R0=R0, a=a_plasma, nfp=nfp, ntheta=ntheta, nzeta=nzeta)
    plasma.net_poloidal_current = 1.0e6
    coil = CoilSurface.circular_torus(R0=R0, a=a_coil, nfp=nfp, ntheta=ntheta, nzeta=nzeta)
    return Regcoil(
        plasma, coil, mpol_potential=mpol, ntor_potential=ntor,
        net_toroidal_current=0.0, symmetry="stellarator_symmetric",
    )


def test_nbf_matches_fourier_mode_count():
    prob = _small_problem(mpol=3, ntor=2)
    # legacy regcoil_init_Fourier_modes(mpol, ntor, ..., include_00=.false.)
    assert prob.nbf == 3 * (2 * 2 + 1) + 2
    assert prob.xm_potential.shape == (prob.nbf,)
    assert prob.xn_potential.shape == (prob.nbf,)
    assert prob.matrix_B.shape == (prob.nbf, prob.nbf)
    assert prob.matrix_K.shape == (prob.nbf, prob.nbf)


def test_symmetry_both_is_twice_stellarator_symmetric():
    sin_only = _small_problem()
    both = Regcoil(
        sin_only.plasma, sin_only.coil, mpol_potential=3, ntor_potential=2,
        net_toroidal_current=0.0, symmetry="both",
    )
    assert both.nbf == 2 * sin_only.nbf


def test_axisymmetric_solution_is_zero_for_every_symmetry_option():
    """A circular-cross-section plasma and coil (different major/minor
    radius) admit an exactly axisymmetric current potential, independent of
    `symmetry` -- a resolution/basis-independent sanity check."""
    plasma = PlasmaSurface.circular_torus(R0=6.0, a=2.0, nfp=3, ntheta=8, nzeta=8)
    plasma.net_poloidal_current = 1.0e6
    coil = CoilSurface.circular_torus(R0=6.0, a=3.0, nfp=3, ntheta=8, nzeta=8)
    for symmetry in ("stellarator_symmetric", "cos_only", "both"):
        prob = Regcoil(
            plasma, coil, mpol_potential=3, ntor_potential=2,
            net_toroidal_current=0.0, symmetry=symmetry,
        )
        sol = prob.solve(1e-10)
        np.testing.assert_allclose(sol.solution, 0, atol=1e-8)
        assert sol.chi2_B < 1e-10


def test_scan_matches_per_lambda_direct_solves():
    prob = _small_problem()
    lambdas = np.array([0.0, 1e-3, 1.0, 1e6, np.inf])
    scanned = prob.scan(lambdas)
    for lam, sol in zip(lambdas, scanned):
        direct = prob.solve(lam)
        np.testing.assert_allclose(sol.solution, direct.solution, rtol=1e-10, atol=1e-12)
        assert sol.chi2_B == pytest.approx(direct.chi2_B, rel=1e-10)
        assert sol.chi2_K == pytest.approx(direct.chi2_K, rel=1e-10)
        assert sol.max_K == pytest.approx(direct.max_K, rel=1e-10)


def test_regcoil_constructor_logs_expensive_steps(caplog):
    with caplog.at_level("INFO"):
        _small_problem()

    messages = [record.getMessage() for record in caplog.records]
    assert any("Starting build_g_and_h" in message for message in messages)
    assert any("Finished build_g_and_h" in message for message in messages)
    assert any("Starting generalized eigensolve" in message for message in messages)
    assert any("Finished generalized eigensolve" in message for message in messages)


def test_two_problems_different_resolutions_do_not_interfere():
    """No shared state: interleaved solves on differently-sized problems
    must not corrupt each other (mirrors tests/unit/test_kernels.py's
    stateless-kernel check, one layer up)."""
    prob_a = _small_problem(ntheta=8, nzeta=8, mpol=3, ntor=2)
    prob_b = _small_problem(ntheta=12, nzeta=10, mpol=5, ntor=4, R0=4.0, a_plasma=1.0, a_coil=1.8, nfp=2)

    sol_a1 = prob_a.solve(1e-3)
    sol_b = prob_b.solve(1e-3)
    sol_a2 = prob_a.solve(1e-3)

    assert prob_a.nbf != prob_b.nbf
    np.testing.assert_array_equal(sol_a1.solution, sol_a2.solution)
    assert sol_a1.solution.shape != sol_b.solution.shape


def test_solve_for_target_matches_direct_solve_at_the_target():
    prob = _small_problem(ntheta=12, nzeta=12, mpol=4, ntor=3)
    lo = prob.solve(0.0).chi2_B
    hi = prob.solve(np.inf).chi2_B
    target = 0.5 * (lo + hi)
    sol = prob.solve_for_target("chi2_B", target)
    assert sol.chi2_B == pytest.approx(target, rel=1e-6)


def test_solve_for_target_raises_outside_achievable_range():
    prob = _small_problem()
    hi = prob.solve(np.inf).chi2_B
    with pytest.raises(ValueError, match="not achievable"):
        prob.solve_for_target("chi2_B", hi * 10 + 1)


def test_current_potential_and_current_density_shapes():
    prob = _small_problem(ntheta=6, nzeta=5)
    sol = prob.solve(1e-3)
    coil = prob.coil
    assert sol.current_potential().shape == (coil.ntheta, coil.nzeta)
    assert sol.current_density().shape == (3, coil.ntheta, coil.nzeta)
    np.testing.assert_array_equal(sol.single_valued_current_potential_mn, sol.solution)


def test_regcoil_rejects_mismatched_nfp():
    plasma = PlasmaSurface.circular_torus(R0=6.0, a=2.0, nfp=3, ntheta=8, nzeta=8)
    plasma.net_poloidal_current = 1.0e6
    coil = CoilSurface.circular_torus(R0=6.0, a=3.0, nfp=4, ntheta=8, nzeta=8)
    with pytest.raises(ValueError, match="nfp"):
        Regcoil(
            plasma, coil, mpol_potential=3, ntor_potential=2,
            symmetry="stellarator_symmetric",
        )


def test_regcoil_defaults_net_poloidal_current_from_plasma():
    plasma = PlasmaSurface.circular_torus(R0=6.0, a=2.0, nfp=3, ntheta=8, nzeta=8)
    plasma.net_poloidal_current = 2.5e6
    coil = CoilSurface.circular_torus(R0=6.0, a=3.0, nfp=3, ntheta=8, nzeta=8)

    prob = Regcoil(
        plasma, coil, mpol_potential=3, ntor_potential=2,
    )

    np.testing.assert_allclose(prob.net_poloidal_current, plasma.net_poloidal_current)