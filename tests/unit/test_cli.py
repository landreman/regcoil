"""Unit tests for the single `regcoil` console-script CLI (Phase 10,
ADR-029): `regcoil plot|compare|cut`, each a thin `load() -> plot/cut ->
show/save` wrapper. Replaces `regcoilPlot`/`compareRegcoil`/
`cutCoilsFromRegcoil`/`cut_saddle_coil`.
"""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import numpy as np
import pytest

import regcoil
from regcoil import CoilSurface, PlasmaSurface, Regcoil
from regcoil._cli import main


@pytest.fixture()
def saved_run(tmp_path):
    plasma = PlasmaSurface.circular_torus(R0=6.0, a=2.0, nfp=3, ntheta=8, nzeta=8)
    plasma.net_poloidal_current = 1.0e6
    coil = CoilSurface.circular_torus(R0=6.0, a=3.0, nfp=3, ntheta=8, nzeta=8)
    prob = Regcoil(plasma, coil, mpol_potential=3, ntor_potential=2)
    scan = prob.scan([0.0, 1e-3, 1.0, np.inf])
    path = tmp_path / "run.nc"
    regcoil.save(path, solutions=scan)
    return path


def test_cli_plot_only_atomic_saves_pdf(saved_run, tmp_path):
    out = tmp_path / "pareto.pdf"
    main(["plot", str(saved_run), "--only", "pareto", "--save", str(out), "--no-show"])
    assert out.exists() and out.stat().st_size > 0


def test_cli_plot_only_plot3d_saves_html(saved_run, tmp_path):
    out = tmp_path / "scene.html"
    main(["plot", str(saved_run), "--only", "plot_3d", "--save", str(out), "--no-show"])
    assert out.exists() and out.stat().st_size > 0


def test_cli_plot_dashboard_saves_multiple_files(saved_run, tmp_path):
    out = tmp_path / "dash.pdf"
    main(["plot", str(saved_run), "--save", str(out), "--no-show"])
    written = list(tmp_path.glob("dash_*"))
    assert len(written) > 1
    assert any(p.suffix == ".html" for p in written)
    assert any(p.suffix == ".pdf" for p in written)


def test_cli_plot_field_map_uses_lambda_index(saved_run, tmp_path):
    out = tmp_path / "cp.pdf"
    main([
        "plot", str(saved_run), "--only", "current_potential",
        "--lambda-index", "2", "--save", str(out), "--no-show",
    ])
    assert out.exists()


def test_cli_compare_overlays_runs(saved_run, tmp_path):
    out = tmp_path / "cmp.pdf"
    main(["compare", str(saved_run), str(saved_run), "--save", str(out), "--no-show"])
    assert out.exists() and out.stat().st_size > 0


def test_cli_cut_writes_makegrid_file(saved_run, tmp_path):
    out = tmp_path / "coils.mytest"
    main([
        "cut", str(saved_run), "--coils-per-half-period", "2",
        "--out", str(out), "--no-show",
    ])
    assert out.exists()
    text = out.read_text()
    assert text.startswith("periods 3")
    assert text.strip().endswith("end")


def test_cli_cut_default_output_name(saved_run, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    main(["cut", str(saved_run), "--coils-per-half-period", "1", "--no-show"])
    produced = list(tmp_path.glob("coils.*"))
    assert len(produced) == 1


def test_cli_requires_a_subcommand():
    with pytest.raises(SystemExit):
        main([])
