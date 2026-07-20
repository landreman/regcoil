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


def _moved_along_normal(plasma, separation, theta, zeta):
    """The moved-along-normal points, computed independently of CoilSurface."""
    ev = plasma._evaluate(theta, zeta)
    r, drdtheta, drdzeta = ev["r"], ev["drdtheta"], ev["drdzeta"]
    normal = np.empty_like(r)
    normal[0] = drdzeta[1] * drdtheta[2] - drdtheta[1] * drdzeta[2]
    normal[1] = drdzeta[2] * drdtheta[0] - drdtheta[2] * drdzeta[0]
    normal[2] = drdzeta[0] * drdtheta[1] - drdtheta[0] * drdzeta[1]
    normal /= np.sqrt(np.sum(normal * normal, axis=0))
    return r + separation * normal


def test_from_uniform_offset_reproduces_moved_points_nonaxisymmetric():
    """The core correctness check for the default (`nu`) construction on a
    genuinely non-axisymmetric plasma: with `mpol`/`ntor` at the transform
    grid's Nyquist limit the DFT->IDFT round-trips, so the reconstructed coil
    surface must pass through exactly the moved-along-normal points. This
    fails for any construction that drops the toroidal-angle shift `nu`
    (X = R*cos(zeta) instead of R*cos(zeta + nu)) on a surface whose normal
    has a toroidal component."""
    plasma = PlasmaSurface(
        xm=[0, 1, 1], xn=[0, 0, 3], rmnc=[5.0, 1.0, 0.15], zmns=[0.0, 1.0, 0.1],
        nfp=3, ntheta=32, nzeta=32,
    )
    sep, nt, nz = 0.4, 16, 14
    coil = CoilSurface.from_uniform_offset(
        plasma, separation=sep, ntheta=nt, nzeta=nz, mpol=nt // 2, ntor=nz // 2
    )

    assert coil.standard_toroidal_angle is False
    assert np.max(np.abs(coil.numns)) > 1e-3  # nu genuinely nonzero here

    theta = 2 * np.pi * np.arange(nt) / nt
    zeta = (2 * np.pi / plasma.nfp) * np.arange(nz) / nz
    r_moved = _moved_along_normal(plasma, sep, theta, zeta)
    np.testing.assert_allclose(coil.r[:, :, :nz], r_moved, atol=1e-9)


def test_from_uniform_offset_nu_has_sine_parity_for_stellarator_symmetric_plasma():
    """A stellarator-symmetric plasma yields an offset surface with only the
    sine (`numns`) angle-shift modes; the cosine part (`numnc`) is zero."""
    plasma = PlasmaSurface(
        xm=[0, 1, 1], xn=[0, 0, 3], rmnc=[5.0, 1.0, 0.15], zmns=[0.0, 1.0, 0.1],
        nfp=3, ntheta=32, nzeta=32,
    )
    assert plasma.stellarator_symmetric
    coil = CoilSurface.from_uniform_offset(plasma, separation=0.4, ntheta=16, nzeta=14, mpol=6, ntor=5)
    assert np.all(coil.numnc == 0)
    assert np.any(coil.numns != 0)
    assert coil.stellarator_symmetric


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
