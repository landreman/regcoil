"""Unit tests for the `Surface` / `FourierSurface` object model."""

import numpy as np
import pytest

from regcoil import FourierSurface


def test_evaluate_hand_checked_case():
    """A 2-mode surface (m=0,n=0 and m=1,n=2 with nfp=2), hand-evaluated at
    three (theta, zetal) points; see the derivation in the PR description.
    """
    xm = [0, 1]
    xn = [0, 2]
    rmnc = [5.0, 1.0]
    rmns = [0.0, 0.5]
    zmnc = [0.0, 0.2]
    zmns = [0.0, 0.3]
    surf = FourierSurface(xm, xn, rmnc, zmns, rmns, zmnc, nfp=2, ntheta=1, nzeta=1)

    theta = np.array([0.0, np.pi / 2])
    zetal = np.array([0.0, np.pi / 2])
    result = surf._evaluate(theta, zetal)
    r = result["r"]
    drdtheta = result["drdtheta"]
    drdzeta = result["drdzeta"]

    # (theta=0, zetal=0): R=6, Z=0.2, phi=0 -> (6, 0, 0.2)
    np.testing.assert_allclose(r[:, 0, 0], [6.0, 0.0, 0.2], atol=1e-12)
    # (theta=pi/2, zetal=0): R=5.5, Z=0.3, phi=0 -> (5.5, 0, 0.3)
    np.testing.assert_allclose(r[:, 1, 0], [5.5, 0.0, 0.3], atol=1e-12)
    # (theta=0, zetal=pi/2): R=4, Z=-0.2, phi=pi/2 -> (0, 4, -0.2)
    np.testing.assert_allclose(r[:, 0, 1], [0.0, 4.0, -0.2], atol=1e-12)

    # Derivatives at (theta=0, zetal=0), hand-derived.
    np.testing.assert_allclose(drdtheta[:, 0, 0], [0.5, 0.0, 0.3], atol=1e-12)
    np.testing.assert_allclose(drdzeta[:, 0, 0], [-1.0, 6.0, -0.6], atol=1e-12)


def test_xn_must_be_multiple_of_nfp():
    """xn includes the factor of nfp (VMEC convention); reject inconsistent input."""
    with pytest.raises(ValueError):
        FourierSurface([1], [1], [1.0], [0.0], nfp=3)


def test_stellarator_symmetric_inferred():
    surf = FourierSurface([0, 1], [0, 0], [1.0, 0.5], [0.0, 0.5], nfp=1)
    assert surf.stellarator_symmetric

    surf2 = FourierSurface([0, 1], [0, 0], [1.0, 0.5], [0.0, 0.5], rmns=[0.0, 0.1], nfp=1)
    assert not surf2.stellarator_symmetric

    # A cosine-parity nu mode (numnc) also breaks stellarator symmetry.
    surf3 = FourierSurface([0, 1], [0, 0], [1.0, 0.5], [0.0, 0.5], numnc=[0.0, 0.1], nfp=1)
    assert not surf3.stellarator_symmetric


def test_nu_angle_shift_places_point_at_zeta_plus_nu():
    """With a single nu mode, the physical toroidal angle is phi = zeta + nu,
    so x = R*cos(zeta + nu), y = R*sin(zeta + nu). Hand-check at a point where
    the shift is nonzero, and confirm the analytic derivatives match finite
    differences."""
    # nu = 0.1 * sin(theta) (m=1, n=0), R = 5 + cos(theta), Z = sin(theta).
    xm, xn = [0, 1], [0, 0]
    surf = FourierSurface(
        xm, xn, rmnc=[5.0, 1.0], zmns=[0.0, 1.0], numns=[0.0, 0.1],
        nfp=1, standard_toroidal_angle=False,
    )
    theta = np.array([0.7])
    zetal = np.array([0.4])
    out = surf._evaluate(theta, zetal)
    R = 5.0 + np.cos(0.7)
    nu = 0.1 * np.sin(0.7)
    phi = 0.4 + nu
    np.testing.assert_allclose(out["r"][:, 0, 0], [R * np.cos(phi), R * np.sin(phi), np.sin(0.7)], atol=1e-12)

    # Analytic vs finite-difference derivatives on a small grid.
    th = np.array([0.3, 1.9, 4.0])
    ze = np.array([0.2, 2.1, 5.0])
    h = 1e-6
    base = surf._evaluate(th, ze)
    fd_dth = (surf._evaluate(th + h, ze)["r"] - surf._evaluate(th - h, ze)["r"]) / (2 * h)
    fd_dze = (surf._evaluate(th, ze + h)["r"] - surf._evaluate(th, ze - h)["r"]) / (2 * h)
    np.testing.assert_allclose(base["drdtheta"], fd_dth, atol=1e-7)
    np.testing.assert_allclose(base["drdzeta"], fd_dze, atol=1e-7)


def test_nu_zero_reduces_to_standard_surface():
    """A FourierSurface with no nu modes evaluates identically whether or not
    the (unused) nu branch would run -- i.e. nu defaults to the phi=zeta fast
    path and standard_toroidal_angle is True."""
    surf = FourierSurface([0, 1, 1], [0, 0, 2], [5.0, 1.0, 0.1], [0.0, 1.0, 0.1], nfp=2)
    assert surf.standard_toroidal_angle is True
    assert np.all(surf.numns == 0) and np.all(surf.numnc == 0)


def test_standard_toroidal_angle_defaults_true():
    """Every plain Fourier-coefficient construction represents its `zeta` as
    the standard toroidal angle; only `CoilSurface.from_uniform_offset`'s
    normal-offset method (default) sets this False."""
    surf = FourierSurface([0, 1], [0, 0], [1.0, 0.5], [0.0, 0.5], nfp=1)
    assert surf.standard_toroidal_angle is True

    surf2 = FourierSurface([0, 1], [0, 0], [1.0, 0.5], [0.0, 0.5], nfp=1, standard_toroidal_angle=False)
    assert surf2.standard_toroidal_angle is False


def test_circular_torus_area_and_shape():
    R0, a, nfp = 5.0, 1.2, 3
    torus = FourierSurface.circular_torus(R0=R0, a=a, nfp=nfp, ntheta=32, nzeta=24)

    assert torus.nfp == nfp
    assert torus.stellarator_symmetric
    np.testing.assert_allclose(torus.area, 4 * np.pi**2 * R0 * a, rtol=1e-12)

    # R, Z on the theta grid at zeta=0 match the standard circular-cross-section formula.
    theta = torus.theta
    R = np.sqrt(torus.r[0, :, 0] ** 2 + torus.r[1, :, 0] ** 2)
    Z = torus.r[2, :, 0]
    np.testing.assert_allclose(R, R0 + a * np.cos(theta), atol=1e-12)
    np.testing.assert_allclose(Z, a * np.sin(theta), atol=1e-12)

    # The legacy volume formula (a 2nd-order-accurate discrete scheme, not a
    # spectral quadrature -- see PR description) converges to the analytic
    # volume 2*pi^2*R0*a^2 as resolution increases; it is not exact at
    # moderate resolution, so check convergence rather than a tight match.
    coarse = FourierSurface.circular_torus(R0=R0, a=a, nfp=nfp, ntheta=16, nzeta=12)
    fine = FourierSurface.circular_torus(R0=R0, a=a, nfp=nfp, ntheta=256, nzeta=192)
    exact = 2 * np.pi**2 * R0 * a**2
    assert abs(fine.volume - exact) < abs(coarse.volume - exact)
    assert abs(fine.volume - exact) < 0.02


def test_cross_section_circular_torus_matches_analytic(nfp=3, R0=5.0, a=1.2):
    """For a plain circular torus, `zeta` already is the physical toroidal
    angle, so cross_section(phi) at any phi should reproduce the exact
    circular cross section (independent of theta-grid resolution)."""
    torus = FourierSurface.circular_torus(R0=R0, a=a, nfp=nfp, ntheta=32, nzeta=16)
    phi = np.array([0.0, 0.3, 1.5, 4.2, 2 * np.pi - 0.01])
    R, Z = torus.cross_section(phi)

    assert R.shape == (len(phi), torus.ntheta)
    assert Z.shape == (len(phi), torus.ntheta)
    # theta values on the grid, in some order, since interpolation at the
    # exact theta grid values (phi is an exact grid zeta) should reproduce
    # them to high precision.
    expected_R = np.sort(R0 + a * np.cos(torus.theta))
    expected_Z = np.sort(a * np.sin(torus.theta))
    for i in range(len(phi)):
        np.testing.assert_allclose(np.sort(R[i]), expected_R, atol=1e-8)
        np.testing.assert_allclose(np.sort(Z[i]), expected_Z, atol=1e-8)


def test_cross_section_default_phi_is_half_field_period():
    torus = FourierSurface.circular_torus(R0=5.0, a=1.0, nfp=4, ntheta=8, nzeta=8)
    R, Z = torus.cross_section()
    assert R.shape == (4, torus.ntheta)


def test_cross_section_nonstandard_toroidal_angle_matches_atan2():
    """A CoilSurface built by `from_uniform_offset(standard_toroidal_angle=False)`
    has a nonzero nu, so a constant-zeta slice is not a constant-phi cross
    section (ADR-025/ADR-026); cross_section must still return the surface's
    actual atan2(y, x) == phi points, not the naive constant-zeta slice."""
    from regcoil import CoilSurface, PlasmaSurface

    plasma = PlasmaSurface.from_wout(
        "equilibria/wout_li383_1.4m.nc", ntheta=12, nzeta=12,
    )
    coil = CoilSurface.from_uniform_offset(
        plasma, separation=0.3, ntheta=12, nzeta=12, mpol=6, ntor=5,
        standard_toroidal_angle=False,
    )
    assert coil.standard_toroidal_angle is False
    assert np.any(coil.numns != 0)

    phi = np.array([0.0, 0.7])
    R, Z = coil.cross_section(phi)
    X = R * np.cos(phi)[:, None]
    Y = R * np.sin(phi)[:, None]
    # Every returned point must actually sit at physical angle phi (up to
    # interpolation error), i.e. atan2(Y, X) == phi.
    recovered_phi = np.arctan2(Y, X)
    for i in range(len(phi)):
        np.testing.assert_allclose(
            np.mod(recovered_phi[i] - phi[i] + np.pi, 2 * np.pi) - np.pi,
            0.0, atol=1e-6,
        )


def test_evaluate_at_matches_grid_evaluate():
    """`FourierSurface.evaluate_at` (arbitrary paired points) must agree with
    `_evaluate` (tensor-product grid) at the same points."""
    surf = FourierSurface.circular_torus(R0=5.0, a=1.0, nfp=3, ntheta=8, nzeta=8)
    theta_grid = surf.theta
    zeta_grid = surf.zeta
    tt, zz = np.meshgrid(theta_grid, zeta_grid, indexing="ij")

    pointwise = surf.evaluate_at(tt.ravel(), zz.ravel())
    grid = surf._evaluate(theta_grid, zeta_grid)

    np.testing.assert_allclose(
        pointwise["r"].reshape(3, surf.ntheta, surf.nzeta), grid["r"], atol=1e-12
    )
    np.testing.assert_allclose(
        pointwise["drdtheta"].reshape(3, surf.ntheta, surf.nzeta), grid["drdtheta"], atol=1e-12
    )
    np.testing.assert_allclose(
        pointwise["drdzeta"].reshape(3, surf.ntheta, surf.nzeta), grid["drdzeta"], atol=1e-12
    )


def test_evaluate_at_derivatives_match_finite_differences_nonstandard_angle():
    """Same check as above, but for a surface with nonzero nu (nonstandard
    toroidal angle), where the derivative formulas have extra terms."""
    from regcoil import CoilSurface, PlasmaSurface

    plasma = PlasmaSurface.circular_torus(R0=5.0, a=1.0, nfp=3, ntheta=8, nzeta=8)
    coil = CoilSurface.from_uniform_offset(
        plasma, separation=0.3, ntheta=8, nzeta=8, mpol=5, ntor=4,
        standard_toroidal_angle=False,
    )
    theta = np.array([0.3, 1.2, 2.5, 4.0])
    zeta = np.array([0.1, 0.4, 0.9, 1.3])
    eps = 1e-6

    base = coil.evaluate_at(theta, zeta)
    dtheta = coil.evaluate_at(theta + eps, zeta)
    dzeta = coil.evaluate_at(theta, zeta + eps)

    np.testing.assert_allclose(
        (dtheta["r"] - base["r"]) / eps, base["drdtheta"], atol=1e-4
    )
    np.testing.assert_allclose(
        (dzeta["r"] - base["r"]) / eps, base["drdzeta"], atol=1e-4
    )
