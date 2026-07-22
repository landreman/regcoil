"""Unit/smoke tests for `regcoil.plot`. Non-interactive
matplotlib backend; checks the atomic-function contract (`ax=`/`fig=` in and
out), that every function works on both a live run and a round-tripped saved
run, and that the saved-run path needs no Fortran kernel / no `eigh`.
"""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pytest

import regcoil
from regcoil import CoilSurface, PlasmaSurface, Regcoil
from regcoil import plot


REPO_ROOT = Path(__file__).resolve().parents[2]
EQUILIBRIA = REPO_ROOT / "equilibria"


def _small_problem(ntheta=8, nzeta=8, mpol=3, ntor=2, nfp=3):
    plasma = PlasmaSurface.circular_torus(R0=6.0, a=2.0, nfp=nfp, ntheta=ntheta, nzeta=nzeta)
    plasma.net_poloidal_current = 1.0e6
    plasma.Bnormal_from_plasma_current = np.random.RandomState(0).randn(ntheta, nzeta) * 1e-3
    coil = CoilSurface.circular_torus(R0=6.0, a=3.0, nfp=nfp, ntheta=ntheta, nzeta=nzeta)
    prob = Regcoil(plasma, coil, mpol_potential=mpol, ntor_potential=ntor)
    return plasma, coil, prob


def _forbid_kernel_calls(monkeypatch):
    import regcoil._core as core
    import scipy.linalg

    def boom(*args, **kwargs):
        raise AssertionError("Fortran kernel / eigh should not be called on a loaded run")

    monkeypatch.setattr(core, "build_g_and_h", boom)
    monkeypatch.setattr(scipy.linalg, "eigh", boom)


@pytest.fixture(autouse=True)
def _close_figures():
    yield
    plt.close("all")


def test_cross_sections_overlay_accepts_and_returns_ax():
    plasma, coil, _ = _small_problem()
    fig, ax = plt.subplots()
    returned = plot.cross_sections_overlay(plasma, ax=ax)
    assert returned is ax

    # No ax= given -> creates its own figure of the documented default size.
    ax2 = plot.cross_sections_overlay(plasma)
    assert ax2.figure.get_size_inches() == pytest.approx(plot.DEFAULT_FIGSIZE)


def test_cross_sections_grid_is_one_subplot_per_phi_red_plasma_blue_coil():
    plasma, coil, _ = _small_problem()
    phi = np.array([0.0, 0.5, 1.0, 1.5])
    fig = plot.cross_sections(plasma, coil, phi=phi)
    assert len(fig.axes) == len(phi)
    for ax in fig.axes:
        colors = {line.get_color() for line in ax.get_lines()}
        assert colors == {"red", "blue"}


def test_pareto_overlays_multiple_scans_on_one_axes():
    _, _, prob = _small_problem()
    scan_a = prob.scan([0.0, 1e-3, 1.0, np.inf])
    scan_b = prob.scan([0.0, 1e-2, 10.0, np.inf])

    fig, ax = plt.subplots()
    returned = plot.pareto([scan_a, scan_b], x="f_K", y="f_B", labels=["a", "b"], ax=ax)
    assert returned is ax
    # Two overlaid lines, one per scan.
    assert len(ax.get_lines()) == 2


def test_lambda_scan_and_atomic_field_maps_return_given_ax():
    _, _, prob = _small_problem()
    scan = prob.scan([0.0, 1e-3, 1.0, np.inf])
    sol = scan[1]

    for fn, kwargs in [
        (plot.lambda_scan, dict()),
        (plot.current_potential, dict(kind="single_valued")),
        (plot.current_potential, dict(kind="total")),
        (plot.current_density, dict()),
        (plot.bnormal, dict(component="plasma_current")),
        (plot.bnormal, dict(component="net_coil")),
        (plot.bnormal, dict(component="total")),
    ]:
        fig, ax = plt.subplots()
        target = sol if fn is not plot.lambda_scan else scan
        returned = fn(target, ax=ax, **kwargs)
        assert returned is ax


def test_current_potential_rejects_unknown_kind():
    _, _, prob = _small_problem()
    sol = prob.solve(1e-3)
    with pytest.raises(ValueError):
        plot.current_potential(sol, kind="bogus")


def test_bnormal_rejects_unknown_component():
    _, _, prob = _small_problem()
    sol = prob.solve(1e-3)
    with pytest.raises(ValueError):
        plot.bnormal(sol, component="bogus")


def test_scan_grids_are_composition_of_atomic_functions():
    _, _, prob = _small_problem()
    scan = prob.scan(np.geomspace(1e-6, 1e6, 10))

    fig1 = plot.current_potential_scan(scan, nmax=4)
    fig2 = plot.bnormal_scan(scan, nmax=4)
    # Each panel contributes a data axes plus its colorbar axes.
    assert len(fig1.axes) == 2 * 4
    assert len(fig2.axes) == 2 * (4 + 2)

    # plot.all() is pure composition too, and works on a lightweight
    # duck-typed container (not just a real regcoil.load() result).
    from types import SimpleNamespace
    plasma, coil, _ = _small_problem()
    data = SimpleNamespace(plasma=plasma, coil=coil, solutions=scan)
    figs = plot.all(data)
    assert len(figs) == 6  # cross-sections grid, pareto, lambda_scan, 2 scan grids, 1 plotly 3d
    assert type(figs[-1]).__module__.startswith("plotly")
    for f in figs[:-1]:
        assert hasattr(f, "savefig")


def test_plot_3d_renders_any_subset():
    plasma, coil, prob = _small_problem()
    sol = prob.solve(1e-3)
    cut = regcoil.cut.cut(sol, coils_per_half_period=2, thickness=0.05)

    fig_plasma_only = plot.plot_3d(plasma=plasma)
    fig_coil_only = plot.plot_3d(winding_surface=coil, winding_surface_style="solid")
    fig_coils_only = plot.plot_3d(coils=cut)
    fig_all = plot.plot_3d(plasma=plasma, winding_surface=coil, coils=cut, winding_surface_style="translucent")

    for fig in (fig_plasma_only, fig_coil_only, fig_coils_only, fig_all):
        assert type(fig).__module__.startswith("plotly")
    assert len(fig_plasma_only.data) > 1
    assert len(fig_all.data) > len(fig_plasma_only.data)
    assert len(fig_all.data) > len(fig_coil_only.data)
    assert len(fig_all.data) > len(fig_coils_only.data)


def test_object_convenience_methods_delegate():
    plasma, coil, prob = _small_problem()
    scan = prob.scan([0.0, 1e-3, np.inf])
    sol = scan[0]

    fig, ax = plt.subplots()
    assert plasma.plot_cross_section(ax=ax) is ax
    assert len(plasma.plot_cross_section(coil).axes) == 4
    assert scan.plot_pareto(ax=ax) is ax
    assert scan.plot_lambda_scan(ax=ax) is ax
    assert sol.plot_current_potential(ax=ax) is ax
    assert sol.plot_current_density(ax=ax) is ax
    assert sol.plot_bnormal(ax=ax) is ax
    assert len(scan.plot_current_potential_scan(nmax=2).axes) == 2 * 2
    assert len(scan.plot_bnormal_scan(nmax=2).axes) == 2 * (2 + 2)


def test_plot_and_cut_need_no_kernel_on_loaded_run(tmp_path, monkeypatch):
    """The ADR-028/029 promise, exercised through the plotting and cutting
    layer: every plot target and cut() work on a saved-then-loaded run
    without calling the Fortran kernel or scipy's generalized eigensolve."""
    plasma, coil, prob = _small_problem()
    scan = prob.scan([0.0, 1e-3, 1.0, np.inf])
    path = tmp_path / "run.nc"
    regcoil.save(path, solutions=scan)

    _forbid_kernel_calls(monkeypatch)

    data = regcoil.load(path)
    sol = data.solutions[1]

    plot.cross_sections(data.plasma, data.coil)
    plot.cross_sections_overlay(data.plasma)
    plot.pareto(data.solutions)
    plot.lambda_scan(data.solutions)
    plot.current_potential(sol, kind="single_valued")
    plot.current_potential(sol, kind="total")
    plot.current_density(sol)
    plot.bnormal(sol, component="plasma_current")
    plot.bnormal(sol, component="net_coil")
    plot.bnormal(sol, component="total")
    plot.current_potential_scan(data.solutions, nmax=4)
    plot.bnormal_scan(data.solutions, nmax=4)
    plot.plot_3d(plasma=data.plasma, winding_surface=data.coil)
    plot.all(data)
    regcoil.cut.cut(sol, coils_per_half_period=2)
