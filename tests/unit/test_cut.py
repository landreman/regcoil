"""Unit tests for `regcoil.cut`: contouring the total
current potential into discrete coils and writing a MAKEGRID file. Replaces
`cutCoilsFromRegcoil` / `cut_saddle_coil`.
"""

from pathlib import Path

import numpy as np
import pytest

import regcoil
from regcoil import CoilSurface, PlasmaSurface, Regcoil
from regcoil.cut import cut


REPO_ROOT = Path(__file__).resolve().parents[2]
EQUILIBRIA = REPO_ROOT / "equilibria"


def _realistic_solution(coils_per_half_period=3, **solve_kwargs):
    """A non-axisymmetric plasma (so the coils are non-trivial), resolved
    enough that contours don't cross theta=0."""
    plasma = PlasmaSurface.from_wout(str(EQUILIBRIA / "wout_d23p4_tm.nc"), ntheta=32, nzeta=32)
    plasma.set_bnormal_from_bnorm_file(str(EQUILIBRIA / "bnorm.d23p4_tm"))
    coil = CoilSurface.from_uniform_offset(
        plasma, separation=0.3, ntheta=32, nzeta=32, mpol=8, ntor=8,
        standard_toroidal_angle=True,
    )
    prob = Regcoil(plasma, coil, mpol_potential=8, ntor_potential=8)
    return prob.solve(**solve_kwargs) if solve_kwargs else prob.solve(1e-13)


def test_cut_produces_expected_number_of_coils():
    sol = _realistic_solution()
    coils_per_half_period = 3
    result = cut(sol, coils_per_half_period=coils_per_half_period)

    expected = 2 * coils_per_half_period * sol.problem.coil.nfp
    assert len(result) == expected
    assert len(result.curves) == expected
    for curve in result.curves:
        assert curve.shape[0] == 3
        # closed: first point repeated at the end.
        np.testing.assert_allclose(curve[:, 0], curve[:, -1])


def test_cut_current_sums_to_net_poloidal_current():
    sol = _realistic_solution()
    result = cut(sol, coils_per_half_period=2)
    total_current = result.current * len(result)
    np.testing.assert_allclose(total_current, sol.problem.net_poloidal_current, rtol=1e-10)


def test_cut_without_thickness_has_no_ribbons():
    sol = _realistic_solution()
    result = cut(sol, coils_per_half_period=2)
    assert result.thickness is None
    with pytest.raises(ValueError):
        result.ribbons()


def test_cut_with_thickness_builds_ribbon_geometry():
    sol = _realistic_solution()
    thickness = 0.05
    result = cut(sol, coils_per_half_period=2, thickness=thickness)
    ribbons = result.ribbons()
    assert len(ribbons) == len(result.curves)
    for corners, curve in zip(ribbons, result.curves):
        assert corners.shape == (4, 3, curve.shape[1])
        # The four corners should be roughly `thickness` apart (the ribbon
        # cross section is a small square of side `thickness`).
        center = corners.mean(axis=0)
        offsets = corners - center[None, :, :]
        dist = np.sqrt(np.sum(offsets**2, axis=1))
        np.testing.assert_allclose(dist, thickness / np.sqrt(2), rtol=0.3)


def test_save_makegrid_writes_valid_file(tmp_path):
    sol = _realistic_solution()
    coils_per_half_period = 2
    result = cut(sol, coils_per_half_period=coils_per_half_period)
    path = tmp_path / "coils.test"
    result.save_makegrid(path)

    lines = path.read_text().splitlines()
    assert lines[0].startswith("periods")
    assert int(lines[0].split()[1]) == sol.problem.coil.nfp
    assert lines[1] == "begin filament"
    assert lines[2] == "mirror NIL"
    assert lines[-1] == "end"

    # One "close the loop" line per coil (5 tokens, last two 0 <coil-id> <name>).
    close_lines = [ln for ln in lines if ln.split()[-1] == "Modular"]
    assert len(close_lines) == len(result.curves)
    for ln in close_lines:
        tokens = ln.split()
        assert float(tokens[3]) == 0.0


def test_cut_needs_no_kernel_on_loaded_run(tmp_path, monkeypatch):
    sol = _realistic_solution()
    path = tmp_path / "run.nc"
    regcoil.save(path, solutions=sol.problem.scan([sol.lam]))

    import regcoil._core as core
    import scipy.linalg

    def boom(*args, **kwargs):
        raise AssertionError("Fortran kernel / eigh should not be called on a loaded run")

    monkeypatch.setattr(core, "build_g_and_h", boom)
    monkeypatch.setattr(scipy.linalg, "eigh", boom)

    data = regcoil.load(path)
    result = cut(data.solutions[0], coils_per_half_period=2, thickness=0.02)
    assert len(result) == 2 * 2 * data.coil.nfp


def test_cut_rejects_theta_wraparound_with_clear_error(monkeypatch):
    """If contour tracing finds more than one path at a level (the coil
    contour crosses theta=0), cut() should raise rather than silently pick
    an arbitrary piece."""
    import contourpy

    sol = _realistic_solution()

    class _FakeGenerator:
        def lines(self, level):
            return [np.zeros((3, 2)), np.zeros((3, 2))]  # two paths instead of one

    monkeypatch.setattr(contourpy, "contour_generator", lambda **kwargs: _FakeGenerator())
    with pytest.raises(RuntimeError, match="theta_shift"):
        cut(sol, coils_per_half_period=2)
