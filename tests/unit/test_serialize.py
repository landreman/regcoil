"""Unit tests for `regcoil.save`/`regcoil.load`: object
serialization to NetCDF-4 via `h5netcdf`. No golden legacy comparison here
(the legacy Fortran output format is not preserved, by design) -- these
tests check the round trip is exact and that a loaded run needs no Fortran
kernel / no `eigh` / no DGEMM to reproduce every plot-target quantity.
"""

from pathlib import Path

import numpy as np
import pytest

import regcoil
from regcoil import CoilSurface, PlasmaSurface, Regcoil, Solution, SolutionScan


REPO_ROOT = Path(__file__).resolve().parents[2]
EQUILIBRIA = REPO_ROOT / "equilibria"


def _small_problem(ntheta=8, nzeta=8, mpol=3, ntor=2, nfp=3, standard_toroidal_angle=True):
    plasma = PlasmaSurface.circular_torus(R0=6.0, a=2.0, nfp=nfp, ntheta=ntheta, nzeta=nzeta)
    plasma.net_poloidal_current = 1.0e6
    plasma.curpol = 1.7
    plasma.Bnormal_from_plasma_current = np.random.RandomState(0).randn(ntheta, nzeta) * 1e-3
    plasma.modB = np.ones((ntheta, nzeta))
    if standard_toroidal_angle:
        coil = CoilSurface.circular_torus(R0=6.0, a=3.0, nfp=nfp, ntheta=ntheta, nzeta=nzeta)
    else:
        coil = CoilSurface.from_uniform_offset(
            plasma, separation=0.5, ntheta=ntheta, nzeta=nzeta, mpol=4, ntor=3,
            standard_toroidal_angle=False,
        )
    prob = Regcoil(
        plasma, coil, mpol_potential=mpol, ntor_potential=ntor,
        net_toroidal_current=0.0, stellarator_symmetric=True,
    )
    return plasma, coil, prob


def _forbid_kernel_calls(monkeypatch):
    """Make the Fortran kernel and the eigendecomposition raise if called,
    so a test using this fixture proves the code path under test needs
    neither (the ADR-028 promise for a loaded run)."""
    import regcoil._core as core
    import scipy.linalg

    def boom(*args, **kwargs):
        raise AssertionError("Fortran kernel / eigh should not be called on a loaded run")

    monkeypatch.setattr(core, "build_g_and_h", boom)
    monkeypatch.setattr(scipy.linalg, "eigh", boom)


def test_round_trip_all_four_object_kinds(tmp_path):
    plasma, coil, prob = _small_problem()
    lambdas = np.array([0.0, 1e-3, 1.0, 1e6, np.inf])
    scan = prob.scan(lambdas)

    path = tmp_path / "run.nc"
    regcoil.save(path, solutions=scan)
    data = regcoil.load(path)

    assert data.plasma is not None and data.coil is not None
    assert data.problem is not None and data.solutions is not None
    assert isinstance(data.solutions, SolutionScan)

    # Surfaces: full Fourier source of truth.
    for name in ("xm", "xn", "rmnc", "rmns", "zmnc", "zmns", "numns", "numnc"):
        np.testing.assert_allclose(getattr(data.plasma, name), getattr(plasma, name))
        np.testing.assert_allclose(getattr(data.coil, name), getattr(coil, name))
    assert data.plasma.nfp == plasma.nfp
    assert data.plasma.stellarator_symmetric == plasma.stellarator_symmetric
    assert data.plasma.standard_toroidal_angle == plasma.standard_toroidal_angle
    np.testing.assert_allclose(data.plasma.net_poloidal_current, plasma.net_poloidal_current)
    np.testing.assert_allclose(data.plasma.curpol, plasma.curpol)
    np.testing.assert_allclose(
        data.plasma.Bnormal_from_plasma_current, plasma.Bnormal_from_plasma_current
    )
    np.testing.assert_allclose(data.plasma.modB, plasma.modB)

    # Problem: scalar params + the one stored operator-derived grid.
    assert data.problem.mpol_potential == prob.mpol_potential
    assert data.problem.ntor_potential == prob.ntor_potential
    np.testing.assert_allclose(data.problem.net_poloidal_current, prob.net_poloidal_current)
    np.testing.assert_allclose(data.problem.net_toroidal_current, prob.net_toroidal_current)
    assert data.problem.stellarator_symmetric == prob.stellarator_symmetric
    np.testing.assert_allclose(
        data.problem.Bnormal_from_net_coil_currents, prob.Bnormal_from_net_coil_currents
    )

    # Every Solution field, for every lambda.
    assert len(data.solutions) == len(scan)
    for loaded, live in zip(data.solutions, scan):
        np.testing.assert_allclose(loaded.lam, live.lam)
        np.testing.assert_allclose(loaded.solution, live.solution)
        np.testing.assert_allclose(loaded.f_B, live.f_B)
        np.testing.assert_allclose(loaded.f_K, live.f_K)
        np.testing.assert_allclose(loaded.max_K, live.max_K)
        np.testing.assert_allclose(loaded.rms_K, live.rms_K)
        np.testing.assert_allclose(loaded.max_Bnormal, live.max_Bnormal)
        np.testing.assert_allclose(loaded.max_Bnormal_over_B, live.max_Bnormal_over_B)
        np.testing.assert_allclose(loaded.avg_Bnormal_over_B, live.avg_Bnormal_over_B)
        np.testing.assert_allclose(loaded.Bnormal_total, live.Bnormal_total)
        np.testing.assert_allclose(loaded.current_potential(), live.current_potential())
        np.testing.assert_allclose(loaded.current_density(), live.current_density())

    np.testing.assert_allclose(data.solutions.lam, scan.lam)
    np.testing.assert_allclose(data.solutions.f_B, scan.f_B)
    np.testing.assert_allclose(data.solutions.f_K, scan.f_K)
    np.testing.assert_allclose(data.solutions.max_K, scan.max_K)
    np.testing.assert_allclose(data.solutions.rms_K, scan.rms_K)
    np.testing.assert_allclose(data.solutions.max_Bnormal, scan.max_Bnormal)
    np.testing.assert_allclose(data.solutions.max_Bnormal_over_B, scan.max_Bnormal_over_B)
    np.testing.assert_allclose(data.solutions.avg_Bnormal_over_B, scan.avg_Bnormal_over_B)


def test_round_trip_needs_no_kernel_no_eigh(tmp_path, monkeypatch):
    """The central ADR-028 promise: plotting a loaded run touches neither
    the Fortran kernel nor scipy's generalized eigensolve."""
    plasma, coil, prob = _small_problem()
    scan = prob.scan(np.array([0.0, 1e-3, np.inf]))
    path = tmp_path / "run.nc"
    regcoil.save(path, solutions=scan)

    _forbid_kernel_calls(monkeypatch)

    data = regcoil.load(path)
    # Surface grids (area/normal/r) -- pure gemm from Fourier coefficients.
    assert data.plasma.area > 0
    assert data.coil.area > 0
    # Every plot-target quantity.
    _ = data.problem.Bnormal_from_net_coil_currents
    for sol in data.solutions:
        _ = sol.Bnormal_total
        _ = sol.current_potential()
        _ = sol.current_density()
    _ = data.solutions.f_B, data.solutions.f_K, data.solutions.max_K, data.solutions.max_Bnormal
    _ = data.solutions.max_Bnormal_over_B, data.solutions.avg_Bnormal_over_B


def test_non_standard_toroidal_angle_coil_round_trips(tmp_path):
    """A CoilSurface with nonzero nu (toroidal-angle-shift) modes round-trips
    only because they are part of the stored source of truth (ADR-026)."""
    plasma = PlasmaSurface.from_wout(str(EQUILIBRIA / "wout_li383_1.4m.nc"), ntheta=8, nzeta=8)
    coil = CoilSurface.from_uniform_offset(
        plasma, separation=0.3, ntheta=8, nzeta=8, mpol=4, ntor=3, standard_toroidal_angle=False
    )
    assert np.any(coil.numns != 0)

    path = tmp_path / "coil.nc"
    coil.save(path)
    loaded = CoilSurface.load(path)
    assert loaded.standard_toroidal_angle is False
    np.testing.assert_allclose(loaded.numns, coil.numns)
    np.testing.assert_allclose(loaded.numnc, coil.numnc)
    np.testing.assert_allclose(loaded.r, coil.r)


def test_single_object_save_load(tmp_path):
    plasma, coil, prob = _small_problem()
    sol = prob.solve(1e-3)

    plasma.save(tmp_path / "plasma.nc")
    p2 = PlasmaSurface.load(tmp_path / "plasma.nc")
    np.testing.assert_allclose(p2.rmnc, plasma.rmnc)
    assert regcoil.load(tmp_path / "plasma.nc").coil is None

    coil.save(tmp_path / "coil.nc")
    c2 = CoilSurface.load(tmp_path / "coil.nc")
    np.testing.assert_allclose(c2.rmnc, coil.rmnc)

    prob.save(tmp_path / "prob.nc")
    p3 = Regcoil.load(tmp_path / "prob.nc")
    assert p3.g is None  # cheap assembly only
    np.testing.assert_allclose(p3.Bnormal_from_net_coil_currents, prob.Bnormal_from_net_coil_currents)
    # A saved-then-loaded problem still supports solving a new lambda.
    np.testing.assert_allclose(p3.solve(1e-3).f_B, sol.f_B, rtol=1e-10)

    sol.save(tmp_path / "sol.nc")
    s2 = Solution.load(tmp_path / "sol.nc")
    np.testing.assert_allclose(s2.solution, sol.solution)
    np.testing.assert_allclose(s2.current_potential(), sol.current_potential())


def test_per_class_load_errors_if_group_missing(tmp_path):
    plasma, coil, prob = _small_problem()
    plasma.save(tmp_path / "plasma_only.nc")

    with pytest.raises(ValueError):
        CoilSurface.load(tmp_path / "plasma_only.nc")
    with pytest.raises(ValueError):
        Regcoil.load(tmp_path / "plasma_only.nc")
    with pytest.raises(ValueError):
        Solution.load(tmp_path / "plasma_only.nc")


def test_solution_load_rejects_multi_lambda_file(tmp_path):
    plasma, coil, prob = _small_problem()
    scan = prob.scan([0.0, 1e-3, 1.0])
    regcoil.save(tmp_path / "scan.nc", solutions=scan)

    with pytest.raises(ValueError):
        Solution.load(tmp_path / "scan.nc")


def test_save_requires_at_least_one_object(tmp_path):
    with pytest.raises(ValueError):
        regcoil.save(tmp_path / "empty.nc")


def test_save_deduplicates_problem_for_a_scan(tmp_path):
    """`solutions=scan` transitively writes one /problem, /plasma, /coil --
    not one per lambda."""
    plasma, coil, prob = _small_problem()
    scan = prob.scan(np.linspace(0.0, 1.0, 5))
    path = tmp_path / "scan.nc"
    regcoil.save(path, solutions=scan)

    import h5netcdf
    with h5netcdf.File(path, "r") as f:
        assert set(f.groups) == {"plasma", "coil", "problem", "solutions"}
        assert f.groups["solutions"].dimensions["lambda"].size == 5
