"""Unit tests for `Regcoil`/`Solution`: object-model behavior that
doesn't need a golden legacy comparison (see tests/regression/ for those).
"""

from pathlib import Path

import numpy as np
import pytest

from regcoil import CoilSurface, PlasmaSurface, Regcoil
from regcoil.regcoil import _build_matrix_K, _flatten_grid


REPO_ROOT = Path(__file__).resolve().parents[2]
EQUILIBRIA = REPO_ROOT / "equilibria"


def _small_problem(ntheta=9, nzeta=8, mpol=3, ntor=2, R0=6.0, a_plasma=2.0, a_coil=3.0, nfp=3):
    plasma = PlasmaSurface.circular_torus(R0=R0, a=a_plasma, nfp=nfp, ntheta=ntheta, nzeta=nzeta)
    plasma.net_poloidal_current = 1.0e6
    plasma.modB = np.ones((ntheta, nzeta))
    coil = CoilSurface.circular_torus(R0=R0, a=a_coil, nfp=nfp, ntheta=ntheta, nzeta=nzeta)
    return Regcoil(
        plasma, coil, mpol_potential=mpol, ntor_potential=ntor,
        net_toroidal_current=0.0, stellarator_symmetric=True,
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
        net_toroidal_current=0.0, stellarator_symmetric=False,
    )
    assert both.nbf == 2 * sin_only.nbf


def test_blocked_matrix_K_matches_direct_contraction():
    rng = np.random.default_rng(1)
    f_all = rng.standard_normal((5, 3, 11))
    norm_normal = rng.uniform(0.5, 2.0, size=11)
    scale = 0.37

    expected = scale * np.einsum(
        "mcg,ncg,g->mn",
        f_all,
        f_all,
        1.0 / norm_normal,
    )
    for grid_block in (1, 4, 100):
        actual = _build_matrix_K(f_all, norm_normal, scale, grid_block=grid_block)
        np.testing.assert_allclose(actual, expected, rtol=2e-15, atol=2e-15)


def test_axisymmetric_solution_vanishes():
    """Circular plasma/coil tori should produce a zero single-valued current
    potential and negligible B_normal at every lambda, independent of
    `stellarator_symmetric`.

    This is a resolution-independent sanity check for the Fourier basis and
    matrix assembly.
    """
    plasma = PlasmaSurface.circular_torus(R0=6.0, a=2.0, nfp=3, ntheta=50, nzeta=48)
    plasma.net_poloidal_current = 1.0e6
    plasma.modB = np.ones((50, 48))
    coil = CoilSurface.circular_torus(R0=6.5, a=4.0, nfp=3, ntheta=54, nzeta=52)
    lambdas = np.geomspace(1e-15, 1e100, num=3)

    for stellarator_symmetric, expected_nbf in ((True, 97), (False, 194)):
        prob = Regcoil(
            plasma, coil, mpol_potential=6, ntor_potential=7,
            net_toroidal_current=0.0,
            stellarator_symmetric=stellarator_symmetric,
        )
        assert prob.nbf == expected_nbf

        for sol in prob.scan(lambdas):
            assert abs(sol.f_B) < 1e-10
            assert abs(sol.max_Bnormal) < 1e-10
            np.testing.assert_allclose(sol.solution, 0, atol=3e-10)


def test_scan_matches_per_lambda_direct_solves():
    prob = _small_problem()
    lambdas = np.array([0.0, 1e-3, 1.0, 1e6, np.inf])
    scanned = prob.scan(lambdas)
    for lam, sol in zip(lambdas, scanned):
        direct = prob.solve(lam)
        np.testing.assert_allclose(sol.solution, direct.solution, rtol=1e-10, atol=1e-12)
        assert sol.f_B == pytest.approx(direct.f_B, rel=1e-10)
        assert sol.f_K == pytest.approx(direct.f_K, rel=1e-10)
        assert sol.max_K == pytest.approx(direct.max_K, rel=1e-10)


def test_bnormal_over_B_diagnostics():
    """`max_Bnormal_over_B` / `avg_Bnormal_over_B` are max and area-mean of
    |B_n / |B|| on the plasma surface."""
    plasma = PlasmaSurface.from_wout(
        str(EQUILIBRIA / "wout_d23p4_tm.nc"), ntheta=17, nzeta=16,
    )
    coil = CoilSurface.from_uniform_offset(
        plasma, separation=0.5, ntheta=19, nzeta=18, standard_toroidal_angle=True,
    )
    prob = Regcoil(plasma, coil, mpol_potential=4, ntor_potential=4)
    sol = prob.solve(1e-3)

    ratio = np.abs(sol.Bnormal_total / plasma.modB)
    assert sol.max_Bnormal_over_B == pytest.approx(float(np.max(ratio)), rel=1e-14)
    expected_avg = float(
        np.sum(ratio * plasma.norm_normal) / np.sum(plasma.norm_normal)
    )
    assert sol.avg_Bnormal_over_B == pytest.approx(expected_avg, rel=1e-14)
    assert sol.max_Bnormal_over_B >= sol.avg_Bnormal_over_B >= 0.0


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
    prob_a = _small_problem(ntheta=8, nzeta=9, mpol=3, ntor=2)
    prob_b = _small_problem(ntheta=12, nzeta=10, mpol=5, ntor=4, R0=4.0, a_plasma=1.0, a_coil=1.8, nfp=2)

    sol_a1 = prob_a.solve(1e-3)
    sol_b = prob_b.solve(1e-3)
    sol_a2 = prob_a.solve(1e-3)

    assert prob_a.nbf != prob_b.nbf
    np.testing.assert_array_equal(sol_a1.solution, sol_a2.solution)
    assert sol_a1.solution.shape != sol_b.solution.shape


def test_solve_for_target_matches_direct_solve_at_the_target():
    prob = _small_problem(ntheta=13, nzeta=12, mpol=4, ntor=3)
    lo = prob.solve(0.0).f_B
    hi = prob.solve(np.inf).f_B
    target = 0.5 * (lo + hi)
    sol = prob.solve_for_target("f_B", target)
    assert sol.f_B == pytest.approx(target, rel=1e-6)


def test_solve_for_target_raises_outside_achievable_range():
    prob = _small_problem()
    hi = prob.solve(np.inf).f_B
    with pytest.raises(ValueError, match="not achievable"):
        prob.solve_for_target("f_B", hi * 10 + 1)


def test_solve_for_target_rejects_unachievable_max_K_targets():
    plasma = PlasmaSurface.from_wout(
        str(EQUILIBRIA / "wout_d23p4_tm.nc"),
        ntheta=65,
        nzeta=64,
    )
    coil = CoilSurface.from_uniform_offset(
        plasma, separation=0.5, ntheta=67, nzeta=66, standard_toroidal_angle=True
    )
    prob = Regcoil(
        plasma, coil, mpol_potential=12, ntor_potential=12,
    )

    with pytest.raises(ValueError, match="not achievable"):
        prob.solve_for_target("max_K", 1.0e6)

    with pytest.raises(ValueError, match="not achievable"):
        prob.solve_for_target("max_K", 2.0e8)


def test_current_potential_and_current_density_shapes():
    prob = _small_problem(ntheta=6, nzeta=5)
    sol = prob.solve(1e-3)
    coil = prob.coil
    assert sol.current_potential().shape == (coil.ntheta, coil.nzeta)
    assert sol.current_density().shape == (3, coil.ntheta, coil.nzeta)
    np.testing.assert_array_equal(sol.single_valued_current_potential_mn, sol.solution)


def _direct_magnetic_field_reference(
    solution, points, *, theta_stride=1, zeta_stride=1
):
    """Temporary-heavy expression retained here as an independent oracle."""
    coil = solution.problem.coil
    K_one_period = solution.current_density()
    source_blocks = []
    current_blocks = []
    weight_blocks = []
    for period in range(coil.nfp):
        start = period * coil.nzeta
        stop = start + coil.nzeta
        source_blocks.append(_flatten_grid(coil.r[:, :, start:stop]).T)

        angle = 2 * np.pi * period / coil.nfp
        c = np.cos(angle)
        s = np.sin(angle)
        K_period = np.empty_like(K_one_period)
        K_period[0] = c * K_one_period[0] - s * K_one_period[1]
        K_period[1] = s * K_one_period[0] + c * K_one_period[1]
        K_period[2] = K_one_period[2]
        current_blocks.append(_flatten_grid(K_period).T)
        weight_blocks.append(
            _flatten_grid(coil.norm_normal) * coil.dtheta * coil.dzeta
        )

    full_grid_shape = (coil.nfp * coil.nzeta, coil.ntheta, 3)
    source_points = np.concatenate(source_blocks).reshape(full_grid_shape)[
        ::zeta_stride, ::theta_stride
    ].reshape(-1, 3)
    surface_current = np.concatenate(current_blocks).reshape(full_grid_shape)[
        ::zeta_stride, ::theta_stride
    ].reshape(-1, 3)
    weights = np.concatenate(weight_blocks).reshape(full_grid_shape[:-1])[
        ::zeta_stride, ::theta_stride
    ].reshape(-1) * theta_stride * zeta_stride
    points = np.asarray(points, dtype=float)
    original_shape = points.shape
    points = points.reshape(-1, 3)
    displacement = points[:, None, :] - source_points[None, :, :]
    distance_squared = np.einsum(
        "ijk,ijk->ij", displacement, displacement
    )
    result = 1.0e-7 * np.einsum(
        "ijk,ij,j->ik",
        np.cross(surface_current[None, :, :], displacement),
        distance_squared ** -1.5,
        weights,
        optimize=True,
    )
    return result.reshape(original_shape)


def test_prepared_magnetic_field_matches_direct_quadrature():
    sol = _small_problem(ntheta=6, nzeta=5).solve(1e-3)
    points = np.array(
        [
            [[6.0, 0.0, 0.0], [6.5, 0.2, -0.1]],
            [[5.5, -0.3, 0.2], [7.0, 0.1, 0.3]],
        ]
    )
    expected = _direct_magnetic_field_reference(sol, points)

    evaluator = sol.prepare_magnetic_field()
    np.testing.assert_allclose(
        evaluator(points, chunk_size=1), expected, rtol=2e-13, atol=2e-13
    )
    np.testing.assert_allclose(
        evaluator(points, chunk_size=3), expected, rtol=2e-13, atol=2e-13
    )
    np.testing.assert_allclose(
        sol.magnetic_field(points), expected, rtol=2e-13, atol=2e-13
    )


def test_magnetic_field_evaluator_is_cached():
    sol = _small_problem(ntheta=6, nzeta=5).solve(1e-3)
    evaluator = sol.prepare_magnetic_field()
    assert sol.prepare_magnetic_field() is evaluator
    sol.magnetic_field([6.0, 0.0, 0.0])
    assert sol.prepare_magnetic_field() is evaluator


def test_single_point_magnetic_field_matches_batched_and_direct_quadrature():
    sol = _small_problem(ntheta=8, nzeta=6, nfp=2).solve(1e-3)
    evaluator = sol.prepare_magnetic_field(theta_stride=2, zeta_stride=2)
    point = np.array([6.5, 0.2, -0.1])
    expected = _direct_magnetic_field_reference(
        sol, point, theta_stride=2, zeta_stride=2
    )

    single = evaluator(point)
    single_from_list = evaluator(point.tolist())
    batched = evaluator(point[None, :])

    assert single.shape == (3,)
    assert single_from_list.shape == (3,)
    assert batched.shape == (1, 3)
    assert evaluator._source_points_c.flags.c_contiguous
    assert evaluator._weighted_surface_current_c.flags.c_contiguous
    np.testing.assert_allclose(single, expected, rtol=2e-13, atol=2e-13)
    np.testing.assert_allclose(single_from_list, expected, rtol=2e-13, atol=2e-13)
    np.testing.assert_allclose(batched[0], expected, rtol=2e-13, atol=2e-13)


def test_reduced_magnetic_field_quadrature_matches_direct_sum():
    sol = _small_problem(ntheta=8, nzeta=6, nfp=2).solve(1e-3)
    points = np.array([[6.0, 0.0, 0.0], [6.5, 0.2, -0.1]])
    expected = _direct_magnetic_field_reference(
        sol, points, theta_stride=2, zeta_stride=2
    )

    evaluator = sol.prepare_magnetic_field(
        theta_stride=2, zeta_stride=2
    )
    assert evaluator.theta_stride == 2
    assert evaluator.zeta_stride == 2
    assert evaluator.quadrature_shape == (6, 4)
    assert evaluator.nsource == 24
    assert evaluator.source_points.shape == (24, 3)
    np.testing.assert_allclose(
        evaluator(points), expected, rtol=2e-13, atol=2e-13
    )
    np.testing.assert_allclose(
        sol.magnetic_field(points, theta_stride=2, zeta_stride=2),
        expected,
        rtol=2e-13,
        atol=2e-13,
    )
    assert sol.prepare_magnetic_field(
        theta_stride=2, zeta_stride=2
    ) is evaluator
    assert sol.prepare_magnetic_field() is not evaluator


@pytest.mark.parametrize(
    "kwargs, error, message",
    [
        ({"theta_stride": 0}, ValueError, "positive"),
        ({"theta_stride": 3}, ValueError, "divide"),
        ({"zeta_stride": 5}, ValueError, "divide"),
        ({"zeta_stride": 1.5}, TypeError, "integer"),
        ({"theta_stride": True}, TypeError, "integer"),
    ],
)
def test_magnetic_field_rejects_invalid_quadrature_strides(
    kwargs, error, message
):
    sol = _small_problem(ntheta=8, nzeta=6, nfp=2).solve(1e-3)
    with pytest.raises(error, match=message):
        sol.prepare_magnetic_field(**kwargs)


@pytest.mark.parametrize(
    "points, message",
    [
        (1.0, "shape"),
        ([1.0, 2.0], "shape"),
        ([np.nan, 0.0, 0.0], "finite"),
    ],
)
def test_magnetic_field_rejects_invalid_points(points, message):
    sol = _small_problem(ntheta=6, nzeta=5).solve(1e-3)
    with pytest.raises(ValueError, match=message):
        sol.magnetic_field(points)


def test_magnetic_field_rejects_nonpositive_chunk_size():
    sol = _small_problem(ntheta=6, nzeta=5).solve(1e-3)
    with pytest.raises(ValueError, match="chunk_size"):
        sol.magnetic_field([6.0, 0.0, 0.0], chunk_size=0)


def test_magnetic_field_rejects_surface_quadrature_point():
    sol = _small_problem(ntheta=6, nzeta=5).solve(1e-3)
    point = sol.problem.coil.r[:, 0, 0]
    with pytest.raises(ValueError, match="quadrature point"):
        sol.magnetic_field(point)


def test_regcoil_rejects_mismatched_nfp():
    plasma = PlasmaSurface.circular_torus(R0=6.0, a=2.0, nfp=3, ntheta=9, nzeta=8)
    plasma.net_poloidal_current = 1.0e6
    coil = CoilSurface.circular_torus(R0=6.0, a=3.0, nfp=4, ntheta=11, nzeta=10)
    with pytest.raises(ValueError, match="nfp"):
        Regcoil(
            plasma, coil, mpol_potential=3, ntor_potential=2,
            stellarator_symmetric=True,
        )


def test_regcoil_defaults_net_poloidal_current_from_plasma():
    plasma = PlasmaSurface.circular_torus(R0=6.0, a=2.0, nfp=3, ntheta=8, nzeta=9)
    plasma.net_poloidal_current = 2.5e6
    coil = CoilSurface.circular_torus(R0=6.0, a=3.0, nfp=3, ntheta=10, nzeta=11)

    prob = Regcoil(
        plasma, coil, mpol_potential=3, ntor_potential=2,
    )

    np.testing.assert_allclose(prob.net_poloidal_current, plasma.net_poloidal_current)
