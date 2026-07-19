"""Unit tests for `CoilSurface`.

`from_nescin` is checked against golden r/normal/area/volume values computed
by the legacy compiled Fortran (`regcoil_init_coil_surface`,
geometry_option_coil=3) on the same nescin file -- see tests/unit/_golden.py.
"""

from pathlib import Path

import numpy as np
import pytest

from regcoil import CoilSurface, PlasmaSurface

from ._golden import NORMAL_COIL, R_COIL

REPO_ROOT = Path(__file__).resolve().parents[2]
EQUILIBRIA = REPO_ROOT / "equilibria"

NESCIN_FILE = EQUILIBRIA / "nescin.w7x_standardConfig_separation0.3"

REF_AREA = 1.9033611117049577e02
REF_VOLUME = 4.9883351694634307e01


def test_from_nescin_matches_legacy():
    coil = CoilSurface.from_nescin(str(NESCIN_FILE), nfp=5, ntheta=4, nzeta=3)

    assert coil.mnmax == 163
    np.testing.assert_allclose(coil.area, REF_AREA, rtol=1e-10)
    np.testing.assert_allclose(coil.volume, REF_VOLUME, rtol=1e-10)
    np.testing.assert_allclose(coil.r, R_COIL, atol=1e-9)
    np.testing.assert_allclose(coil.normal, NORMAL_COIL, atol=1e-9)


def test_filter_modes_zeroes_high_harmonics():
    coil = CoilSurface.from_nescin(str(NESCIN_FILE), nfp=5, ntheta=8, nzeta=8)
    assert np.any(np.abs(coil.xm) > 2) or np.any(np.abs(coil.xn) > 2 * coil.nfp)

    coil.filter_modes(mpol_filter=2, ntor_filter=2)
    kept = (np.abs(coil.xm) <= 2) & (np.abs(coil.xn) <= 2 * coil.nfp)
    assert np.all(coil.rmnc[~kept] == 0)
    assert np.all(coil.rmns[~kept] == 0)
    assert np.all(coil.zmnc[~kept] == 0)
    assert np.all(coil.zmns[~kept] == 0)


def test_from_nescin_filter_at_construction():
    filtered = CoilSurface.from_nescin(str(NESCIN_FILE), nfp=5, ntheta=8, nzeta=8, mpol_filter=2, ntor_filter=2)
    unfiltered = CoilSurface.from_nescin(str(NESCIN_FILE), nfp=5, ntheta=8, nzeta=8)
    assert not np.allclose(filtered.rmnc, unfiltered.rmnc)


@pytest.mark.parametrize("standard_toroidal_angle", [False, True])
def test_from_uniform_offset_circular_torus_is_exact(standard_toroidal_angle):
    """`from_uniform_offset` for a circular-cross-section plasma should stay
    circular, with minor radius a + separation, under *either* construction:
    for an axisymmetric plasma the normal has no toroidal component, so the
    default normal-offset method and the standard-toroidal-angle root solve
    agree exactly (see test_kernels.py for the golden-legacy-Fortran and
    non-circular-plasma coverage)."""
    plasma = PlasmaSurface.circular_torus(R0=5.0, a=1.0, nfp=3, ntheta=8, nzeta=8)
    coil = CoilSurface.from_uniform_offset(
        plasma, separation=0.3, ntheta=8, nzeta=8, mpol=2, ntor=1,
        standard_toroidal_angle=standard_toroidal_angle,
    )

    assert coil.standard_toroidal_angle is standard_toroidal_angle
    for m, n, rc, zs in zip(coil.xm, coil.xn, coil.rmnc, coil.zmns):
        if m == 0 and n == 0:
            assert rc == pytest.approx(5.0, abs=1e-9)
        elif m == 1 and n == 0:
            assert rc == pytest.approx(1.3, abs=1e-9)
            assert zs == pytest.approx(1.3, abs=1e-9)
        else:
            assert rc == pytest.approx(0.0, abs=1e-8)
            assert zs == pytest.approx(0.0, abs=1e-8)


def test_from_uniform_offset_default_is_not_standard_toroidal_angle():
    plasma = PlasmaSurface.circular_torus(R0=5.0, a=1.0, nfp=3, ntheta=8, nzeta=8)
    coil = CoilSurface.from_uniform_offset(plasma, separation=0.3, ntheta=8, nzeta=8, mpol=2, ntor=1)
    assert coil.standard_toroidal_angle is False


def test_from_uniform_offset_logs_kernel_timing(caplog):
    plasma = PlasmaSurface.circular_torus(R0=5.0, a=1.0, nfp=3, ntheta=8, nzeta=8)

    with caplog.at_level("INFO"):
        CoilSurface.from_uniform_offset(
            plasma, separation=0.3, ntheta=8, nzeta=8, mpol=2, ntor=1, standard_toroidal_angle=True
        )

    messages = [record.getMessage() for record in caplog.records]
    assert any("Starting uniform offset surface kernel" in message for message in messages)
    assert any("Finished uniform offset surface kernel" in message for message in messages)


def test_from_uniform_offset_default_does_not_log_kernel_timing(caplog):
    """The default (`standard_toroidal_angle=False`) construction is pure
    Python/numpy and never calls the Fortran kernel, so it must not emit the
    kernel's timing log lines."""
    plasma = PlasmaSurface.circular_torus(R0=5.0, a=1.0, nfp=3, ntheta=8, nzeta=8)

    with caplog.at_level("INFO"):
        CoilSurface.from_uniform_offset(plasma, separation=0.3, ntheta=8, nzeta=8, mpol=2, ntor=1)

    messages = [record.getMessage() for record in caplog.records]
    assert not any("uniform offset surface kernel" in message for message in messages)
