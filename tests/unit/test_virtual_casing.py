"""Unit tests for reading B_normal from simsopt virtual-casing data.

simsopt is not a dependency of REGCOIL and is not installed in CI, so the
virtual-casing NetCDF files used here are written by `write_vcasing_file`
below, which reproduces the subset of `simsopt.mhd.VirtualCasing.save`'s output
that REGCOIL reads. The grids follow simsopt's `Surface.get_phi_quadpoints`
exactly (see `simsopt_phi_grid`).

The main test drives the data from the bundled `bnorm.li383_1.4m` file: the
same B_normal is fed in through both the BNORM and the virtual-casing path, and
the two must agree. (That the two formats really do hold the same quantity, in
the same sign convention, is simsopt's
`tests/mhd/test_virtual_casing.py::test_bnorm_benchmark`.)
"""

from pathlib import Path

import numpy as np
import pytest
from scipy.io import netcdf_file

from regcoil import PlasmaSurface, examples

REPO_ROOT = Path(__file__).resolve().parents[2]
EQUILIBRIA = REPO_ROOT / "equilibria"

WOUT_LI383 = EQUILIBRIA / "wout_li383_1.4m.nc"
BNORM_LI383 = EQUILIBRIA / "bnorm.li383_1.4m"


def simsopt_phi_grid(nphi, nfp, use_stellsym):
    """simsopt's toroidal quadrature points, in units where the period is 1.

    `range="half period"` for a stellarator-symmetric calculation (half a field
    period, shifted by half a grid spacing), `range="field period"` otherwise.
    """
    end = (0.5 if use_stellsym else 1.0) / nfp
    phi = np.linspace(0.0, end, nphi, endpoint=False)
    if use_stellsym:
        phi = phi + 0.5 * (phi[1] - phi[0])
    return phi


def write_vcasing_file(path, nfp, trgt_phi, trgt_theta, B_external_normal):
    """Write the variables REGCOIL reads from a simsopt vcasing*.nc file."""
    with netcdf_file(str(path), "w") as f:
        f.createDimension("trgt_ntheta", len(trgt_theta))
        f.createDimension("trgt_nphi", len(trgt_phi))
        for name, value in (
            ("nfp", nfp),
            ("trgt_ntheta", len(trgt_theta)),
            ("trgt_nphi", len(trgt_phi)),
        ):
            var = f.createVariable(name, "i", tuple())
            var.data[()] = value
        f.createVariable("trgt_theta", "d", ("trgt_ntheta",))[:] = trgt_theta
        f.createVariable("trgt_phi", "d", ("trgt_nphi",))[:] = trgt_phi
        f.createVariable("B_external_normal", "d", ("trgt_nphi", "trgt_ntheta"))[
            :, :
        ] = B_external_normal


def bnorm_series(filename, theta, phi, nfp, curpol):
    """Evaluate a BNORM file's Fourier series at *paired* (theta, phi) arrays,
    with theta and phi in radians. This is the same sum as
    `_io.read_bnorm_file`, but written independently and evaluated on simsopt's
    grid convention rather than REGCOIL's.
    """
    modes = np.loadtxt(filename)
    Bn = np.zeros(np.broadcast(theta, phi).shape)
    for m, n, amplitude in modes:
        Bn += amplitude * np.sin(m * theta + n * nfp * phi)
    return Bn * curpol


def li383_vcasing(path, nphi, ntheta, use_stellsym, nfp=3):
    """Write a vcasing file holding the bundled li383 BNORM data, sampled on
    simsopt's target grid. Returns the file path.
    """
    curpol = PlasmaSurface.from_wout(str(WOUT_LI383), ntheta=4, nzeta=4).curpol
    phi = simsopt_phi_grid(nphi, nfp, use_stellsym)
    theta = np.linspace(0.0, 1.0, ntheta, endpoint=False)
    B = bnorm_series(
        BNORM_LI383, 2 * np.pi * theta[None, :], 2 * np.pi * phi[:, None], nfp, curpol
    )
    write_vcasing_file(path, nfp, phi, theta, B)
    return path


# The bundled BNORM file has |m| <= 24 and |n| <= 24, so a 64-point theta grid
# and a 64-point full-field-period phi grid (32 points on the half period that
# stellarator symmetry extends to 64) resolve it exactly. The interpolation is
# then exact too, and the two input paths must agree to roundoff.
NPHI = 32
NTHETA = 64


@pytest.mark.parametrize("use_stellsym", [True, False])
@pytest.mark.parametrize("ntheta,nzeta", [(5, 4), (16, 13), (64, 32)])
def test_matches_bnorm_file(tmp_path, use_stellsym, ntheta, nzeta):
    nphi = NPHI if use_stellsym else 2 * NPHI
    vcasing = li383_vcasing(tmp_path / "vcasing.nc", nphi, NTHETA, use_stellsym)

    reference = PlasmaSurface.from_wout(str(WOUT_LI383), ntheta=ntheta, nzeta=nzeta)
    reference.set_bnormal_from_bnorm_file(str(BNORM_LI383))

    plasma = PlasmaSurface.from_wout(str(WOUT_LI383), ntheta=ntheta, nzeta=nzeta)
    plasma.set_bnormal_from_virtual_casing(str(vcasing))

    assert plasma.Bnormal_from_plasma_current.shape == (ntheta, nzeta)
    np.testing.assert_allclose(
        plasma.Bnormal_from_plasma_current,
        reference.Bnormal_from_plasma_current,
        atol=1e-11,
    )


def test_stellsym_and_non_stellsym_layouts_agree(tmp_path):
    """The same field written on the half-period and the full-period grid must
    interpolate to the same values."""
    half = li383_vcasing(tmp_path / "half.nc", NPHI, NTHETA, use_stellsym=True)
    full = li383_vcasing(tmp_path / "full.nc", 2 * NPHI, NTHETA, use_stellsym=False)

    results = []
    for filename in (half, full):
        plasma = PlasmaSurface.from_wout(str(WOUT_LI383), ntheta=17, nzeta=11)
        plasma.set_bnormal_from_virtual_casing(str(filename))
        results.append(plasma.Bnormal_from_plasma_current)

    np.testing.assert_allclose(results[0], results[1], atol=1e-11)


def test_accepts_virtual_casing_object(tmp_path):
    """A VirtualCasing-like object works as well as a file. simsopt's
    `VirtualCasing.load` leaves scalars as 0-d arrays, so nfp is passed that way
    here."""

    class FakeVirtualCasing:
        pass

    nfp = 3
    phi = simsopt_phi_grid(NPHI, nfp, use_stellsym=True)
    theta = np.linspace(0.0, 1.0, NTHETA, endpoint=False)
    curpol = PlasmaSurface.from_wout(str(WOUT_LI383), ntheta=4, nzeta=4).curpol

    vc = FakeVirtualCasing()
    vc.nfp = np.array(nfp)
    vc.trgt_phi = phi
    vc.trgt_theta = theta
    vc.B_external_normal = bnorm_series(
        BNORM_LI383, 2 * np.pi * theta[None, :], 2 * np.pi * phi[:, None], nfp, curpol
    )

    plasma = PlasmaSurface.from_wout(str(WOUT_LI383), ntheta=9, nzeta=7)
    plasma.set_bnormal_from_virtual_casing(vc)

    reference = PlasmaSurface.from_wout(str(WOUT_LI383), ntheta=9, nzeta=7)
    reference.set_bnormal_from_bnorm_file(str(BNORM_LI383))
    np.testing.assert_allclose(
        plasma.Bnormal_from_plasma_current,
        reference.Bnormal_from_plasma_current,
        atol=1e-11,
    )


def test_interpolation_is_exact_on_the_simsopt_grid(tmp_path):
    """Where the REGCOIL and simsopt grids coincide, the interpolant must
    return the stored values themselves.

    The non-stellsym layout is the one whose phi grid can coincide with
    REGCOIL's: it starts at phi=0 and spans a field period. The half-period
    layout never can, since it is offset by half a grid spacing.
    """
    nfp, nphi, ntheta = 3, 12, 10
    rng = np.random.default_rng(0)
    phi = simsopt_phi_grid(nphi, nfp, use_stellsym=False)
    theta = np.linspace(0.0, 1.0, ntheta, endpoint=False)
    B = rng.standard_normal((nphi, ntheta))
    path = tmp_path / "v.nc"
    write_vcasing_file(path, nfp, phi, theta, B)

    plasma = PlasmaSurface.from_wout(str(WOUT_LI383), ntheta=ntheta, nzeta=nphi)
    plasma.set_bnormal_from_virtual_casing(str(path))

    np.testing.assert_allclose(plasma.Bnormal_from_plasma_current, B.T, atol=1e-12)


def test_nfp_mismatch_raises(tmp_path):
    path = tmp_path / "v.nc"
    nfp = 5  # li383 has nfp=3
    write_vcasing_file(
        path,
        nfp,
        simsopt_phi_grid(8, nfp, use_stellsym=True),
        np.linspace(0.0, 1.0, 8, endpoint=False),
        np.zeros((8, 8)),
    )

    plasma = PlasmaSurface.from_wout(str(WOUT_LI383), ntheta=8, nzeta=8)
    with pytest.raises(ValueError, match="nfp"):
        plasma.set_bnormal_from_virtual_casing(str(path))


def test_bad_source_type_raises():
    plasma = PlasmaSurface.from_wout(str(WOUT_LI383), ntheta=4, nzeta=4)
    with pytest.raises(TypeError, match="B_external_normal"):
        plasma.set_bnormal_from_virtual_casing(object())


@pytest.mark.parametrize("name", ["li383_1.4m", "d23p4_tm"])
def test_matches_bnorm_file_for_bundled_examples(name):
    """The bundled `vcasing*.nc` files (real simsopt output, one per example
    dataset) and the bundled BNORM files describe the same plasma boundaries,
    so B_normal loaded from each should agree.

    Unlike `test_matches_bnorm_file` above, these are two independent physical
    calculations (the legacy BNORM code vs. simsopt's virtual-casing solve),
    not the same series pushed through two code paths, so exact agreement
    isn't expected -- only that they describe the same field to within the
    discretization/method error between the two calculations.
    """
    ds = examples(name)

    reference = PlasmaSurface.from_wout(str(ds.wout), ntheta=32, nzeta=32)
    reference.set_bnormal_from_bnorm_file(str(ds.bnorm))

    plasma = PlasmaSurface.from_wout(str(ds.wout), ntheta=32, nzeta=32)
    plasma.set_bnormal_from_virtual_casing(str(ds.vcasing))

    assert plasma.Bnormal_from_plasma_current.shape == (32, 32)
    rms = np.sqrt(np.mean(reference.Bnormal_from_plasma_current**2))

    # import matplotlib.pyplot as plt
    # plt.figure(figsize=(14.5, 8.1))
    # plt.subplot(2, 2, 1)
    # plt.contourf(
    #     reference.Bnormal_from_plasma_current,
    #     levels=20,
    #     cmap="RdBu_r",
    #     extend="both",
    # )
    # plt.colorbar()
    # plt.title("BNORM")

    # plt.subplot(2, 2, 2)
    # plt.contourf(
    #     plasma.Bnormal_from_plasma_current,
    #     levels=20,
    #     cmap="RdBu_r",
    #     extend="both",
    # )
    # plt.colorbar()
    # plt.title("Virtual Casing")

    # plt.subplot(2, 2, 3)
    # plt.contourf(
    #     plasma.Bnormal_from_plasma_current - reference.Bnormal_from_plasma_current,
    #     levels=20,
    #     cmap="RdBu_r",
    #     extend="both",
    # )
    # plt.colorbar()
    # plt.title("Difference (VC - BNORM)")

    # plt.tight_layout()
    # plt.show()

    np.testing.assert_allclose(
        plasma.Bnormal_from_plasma_current,
        reference.Bnormal_from_plasma_current,
        atol=0.19 * rms,
    )
