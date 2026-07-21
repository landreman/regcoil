"""Unit tests for `regcoil.reparameterize` (ADR-031).

The defining property of a theta reparameterization is that it changes the
*labels* of a surface's points but not the surface, so most of these tests are
invariance checks: the same physical geometry, differently parameterized.
"""

from pathlib import Path

import numpy as np
import pytest

from regcoil import CoilSurface, CurvatureWeighted, PlasmaSurface, UniformArclength
from regcoil.reparameterize import (
    ThetaMap,
    _as_scheme,
    _cumulative_integral,
    _fft_derivative,
    theta_map,
)

REPO_ROOT = Path(__file__).resolve().parents[2]
EQUILIBRIA = REPO_ROOT / "equilibria"

# A markedly eccentric elliptical cross section: the naive angle's arclength
# varies by a factor of ~a/b, so a correct map has real work to do.
R0, A_MAJOR, B_MINOR = 5.0, 1.0, 0.3


def _ellipse_curve(theta, zeta):
    """A torus of elliptical cross section, as a `theta_map` curve callable."""
    R = R0 + A_MAJOR * np.cos(theta)[:, None] + 0 * zeta[None, :]
    Z = B_MINOR * np.sin(theta)[:, None] + 0 * zeta[None, :]
    phi = zeta[None, :] + 0 * theta[:, None]
    return np.stack([R * np.cos(phi), R * np.sin(phi), Z], axis=0)


def _arclengths(theta_old, zeta):
    """Chord lengths between consecutive points of the elliptical curve."""
    R = R0 + A_MAJOR * np.cos(theta_old)
    Z = B_MINOR * np.sin(theta_old)
    pts = np.stack([R * np.cos(zeta), R * np.sin(zeta), Z], axis=0)
    return np.sqrt(np.sum((np.roll(pts, -1, axis=1) - pts) ** 2, axis=0))


# --------------------------------------------------------------------------
# Spectral primitives
# --------------------------------------------------------------------------


def test_fft_derivative_matches_analytic():
    n = 64
    theta = 2 * np.pi * np.arange(n) / n
    # Layout (component, theta, zeta), as `_weight` uses.
    f = np.stack([np.cos(3 * theta) + 0.5 * np.sin(theta), np.exp(np.cos(theta))])[:, :, None]
    d1 = np.stack(
        [-3 * np.sin(3 * theta) + 0.5 * np.cos(theta), -np.sin(theta) * np.exp(np.cos(theta))]
    )[:, :, None]
    d2 = np.stack(
        [
            -9 * np.cos(3 * theta) - 0.5 * np.sin(theta),
            (np.sin(theta) ** 2 - np.cos(theta)) * np.exp(np.cos(theta)),
        ]
    )[:, :, None]

    np.testing.assert_allclose(_fft_derivative(f, 1, axis=1), d1, atol=1e-12)
    np.testing.assert_allclose(_fft_derivative(f, 2, axis=1), d2, atol=1e-11)


def test_cumulative_integral_matches_analytic_and_is_anchored():
    n = 64
    theta = 2 * np.pi * np.arange(n) / n
    w = (2.0 + np.cos(theta) + 0.3 * np.sin(2 * theta))[:, None]
    expected = (2 * theta + np.sin(theta) + 0.3 * (1 - np.cos(2 * theta)) / 2)[:, None]

    result = _cumulative_integral(w)
    np.testing.assert_allclose(result, expected, atol=1e-13)
    # F(0) = 0 must hold *exactly*, not just to roundoff: the periodic spline
    # that inverts the map rejects endpoints that differ at all.
    assert result[0, 0] == 0.0


# --------------------------------------------------------------------------
# The map itself
# --------------------------------------------------------------------------


def test_uniform_arclength_equalizes_arclength_on_an_ellipse():
    ntheta, nzeta = 128, 4
    tmap = theta_map(_ellipse_curve, UniformArclength(), nfp=3, ntheta=ntheta, nzeta=nzeta)

    assert tmap.diagnostics["max_residual"] < 1e-9

    theta_out = 2 * np.pi * np.arange(ntheta) / ntheta
    zeta = (2 * np.pi / 3) * np.arange(nzeta) / nzeta
    theta_old = tmap(theta_out, zeta)

    # In the *old* angle the arclength varies by more than 3x; in the new one
    # the chord lengths agree to the O(h^2) chord-vs-arc error.
    naive = _arclengths(np.broadcast_to(theta_out[:, None], theta_old.shape), zeta)
    assert naive.max() / naive.min() > 3.0
    mapped = _arclengths(theta_old, zeta)
    assert (mapped.max(axis=0) - mapped.min(axis=0)).max() / mapped.mean() < 1e-2


def test_theta_map_is_monotone_periodic_and_anchored():
    ntheta, nzeta = 64, 4
    tmap = theta_map(_ellipse_curve, UniformArclength(), nfp=3, ntheta=ntheta, nzeta=nzeta)
    theta_out = 2 * np.pi * np.arange(ntheta) / ntheta
    zeta = (2 * np.pi / 3) * np.arange(nzeta) / nzeta
    theta_old = tmap(theta_out, zeta)

    # Anchored at theta = 0 (ADR-031 decision 6), strictly increasing, and
    # advancing by exactly 2*pi over one poloidal turn.
    np.testing.assert_allclose(theta_old[0], 0.0, atol=1e-12)
    assert np.all(np.diff(theta_old, axis=0) > 0)
    wrapped = tmap(theta_out + 2 * np.pi, zeta)
    np.testing.assert_allclose(wrapped, theta_old + 2 * np.pi, atol=1e-10)


def test_curvature_weighted_makes_arclength_track_inverse_curvature():
    """`dl/dtheta` proportional to `kappa^-exponent` is the defining property.
    On the ellipse, curvature is highest at the ends of the major axis
    (theta = 0, pi), so those get the *shortest* steps -- the opposite of what
    uniform arclength does."""
    ntheta, nzeta = 256, 2
    scheme = CurvatureWeighted(exponent=1.0, floor=0.0)
    tmap = theta_map(_ellipse_curve, scheme, nfp=3, ntheta=ntheta, nzeta=nzeta)
    theta_out = 2 * np.pi * np.arange(ntheta) / ntheta
    zeta = (2 * np.pi / 3) * np.arange(nzeta) / nzeta
    theta_old = tmap(theta_out, zeta)

    dl = _arclengths(theta_old, zeta)
    # Analytic curvature of (a cos t, b sin t), evaluated at each segment's
    # *midpoint*. Using an endpoint instead leaves an O(h) error that swamps
    # the property under test, since kappa varies by ~10x around this ellipse.
    ahead = np.roll(theta_old, -1, axis=0)
    ahead = ahead + np.where(ahead < theta_old, 2 * np.pi, 0.0)
    t = 0.5 * (theta_old + ahead)
    kappa = (A_MAJOR * B_MINOR) / (
        (A_MAJOR**2 * np.sin(t) ** 2 + B_MINOR**2 * np.cos(t) ** 2) ** 1.5
    )
    product = dl * kappa
    np.testing.assert_allclose(product, product.mean(), rtol=2e-2)

    # And the ordering is genuinely inverted relative to uniform arclength.
    uniform = theta_map(_ellipse_curve, UniformArclength(), nfp=3, ntheta=ntheta, nzeta=nzeta)
    dl_uniform = _arclengths(uniform(theta_out, zeta), zeta)
    assert dl[0, 0] < dl_uniform[0, 0]


def test_curvature_floor_regularizes_a_straight_segment():
    """With `floor=0` a curve with an inflection (kappa -> 0) allots unbounded
    arclength there; the floor must tame it."""
    ntheta = 128

    def curve(theta, zeta):
        # A cross section with a near-straight stretch: kappa passes near zero.
        R = R0 + np.cos(theta)[:, None] + 0.45 * np.cos(2 * theta)[:, None] + 0 * zeta[None, :]
        Z = np.sin(theta)[:, None] + 0.45 * np.sin(2 * theta)[:, None] + 0 * zeta[None, :]
        phi = zeta[None, :] + 0 * theta[:, None]
        return np.stack([R * np.cos(phi), R * np.sin(phi), Z], axis=0)

    unfloored = theta_map(
        curve, CurvatureWeighted(exponent=1.0, floor=0.0), nfp=2, ntheta=ntheta, nzeta=2
    )
    floored = theta_map(
        curve, CurvatureWeighted(exponent=1.0, floor=0.2), nfp=2, ntheta=ntheta, nzeta=2
    )

    theta_out = 2 * np.pi * np.arange(ntheta) / ntheta
    zeta = np.array([0.0])
    spread = lambda tm: np.diff(tm(theta_out, zeta), axis=0).max() / np.diff(  # noqa: E731
        tm(theta_out, zeta), axis=0
    ).min()
    assert spread(floored) < spread(unfloored)


def test_uniform_arclength_is_identity_on_a_circular_cross_section():
    """A circular cross section already has uniform arclength, so the map must
    be (numerically) the identity -- `g - theta` has no Fourier content."""
    ntheta, nzeta = 64, 4

    def circle(theta, zeta):
        R = R0 + np.cos(theta)[:, None] + 0 * zeta[None, :]
        Z = np.sin(theta)[:, None] + 0 * zeta[None, :]
        phi = zeta[None, :] + 0 * theta[:, None]
        return np.stack([R * np.cos(phi), R * np.sin(phi), Z], axis=0)

    tmap = theta_map(circle, UniformArclength(), nfp=3, ntheta=ntheta, nzeta=nzeta)
    theta_out = 2 * np.pi * np.arange(ntheta) / ntheta
    zeta = (2 * np.pi / 3) * np.arange(nzeta) / nzeta

    np.testing.assert_allclose(
        tmap(theta_out, zeta), np.broadcast_to(theta_out[:, None], (ntheta, nzeta)), atol=1e-10
    )
    assert np.abs(tmap.gmns).max() < 1e-10


def test_measure_3d_and_poloidal_agree_when_nu_is_zero():
    """The two measures differ only by the `R^2 nu'^2` term, which vanishes
    when `zeta` is the standard toroidal angle -- as it is for `_ellipse_curve`."""
    kwargs = dict(nfp=3, ntheta=64, nzeta=4)
    a = theta_map(_ellipse_curve, UniformArclength(measure="3d"), **kwargs)
    b = theta_map(_ellipse_curve, UniformArclength(measure="poloidal"), **kwargs)
    np.testing.assert_allclose(a.gmns, b.gmns, atol=1e-12)


def test_stellarator_symmetry_gives_the_map_sine_parity():
    """`g - theta` is odd for a stellarator-symmetric surface (a consequence of
    the theta=0 anchor), so the cosine modes vanish."""
    plasma = PlasmaSurface.from_wout(str(EQUILIBRIA / "wout_d23p4_tm.nc"), ntheta=8, nzeta=8)
    tmap = theta_map(
        lambda th, ze: plasma._evaluate(th, ze)["r"],
        UniformArclength(), nfp=plasma.nfp, ntheta=32, nzeta=16,
        stellarator_symmetric=False,  # do not *impose* the parity; measure it
    )
    assert np.abs(tmap.gmnc).max() < 1e-10
    assert np.abs(tmap.gmns).max() > 1e-3


def test_scheme_shorthand_strings_and_errors():
    assert _as_scheme("uniform_arclength") == UniformArclength()
    assert _as_scheme("curvature") == CurvatureWeighted()
    scheme = CurvatureWeighted(exponent=0.5)
    assert _as_scheme(scheme) is scheme

    with pytest.raises(ValueError, match="Unknown theta reparameterization"):
        _as_scheme("nonsense")
    with pytest.raises(TypeError, match="must be a UniformArclength"):
        _as_scheme(3.0)


def test_theta_map_rejects_a_badly_shaped_curve():
    with pytest.raises(ValueError, match="must return shape"):
        theta_map(
            lambda th, ze: np.zeros((3, len(th), len(ze) + 1)),
            UniformArclength(), nfp=3, ntheta=16, nzeta=4,
        )


def test_theta_map_is_serializable_as_fourier_modes():
    """The map is stored as modes (ADR-031 decision 8), and at the default
    (Nyquist) resolution reproduces its grid values exactly."""
    ntheta, nzeta = 32, 8
    tmap = theta_map(_ellipse_curve, UniformArclength(), nfp=3, ntheta=ntheta, nzeta=nzeta)
    assert isinstance(tmap, ThetaMap)
    assert tmap.scheme == UniformArclength()

    rebuilt = ThetaMap(tmap.xm, tmap.xn, tmap.gmns, tmap.gmnc, nfp=tmap.nfp)
    theta_out = 2 * np.pi * np.arange(ntheta) / ntheta
    zeta = (2 * np.pi / 3) * np.arange(nzeta) / nzeta
    np.testing.assert_allclose(rebuilt(theta_out, zeta), tmap(theta_out, zeta), atol=1e-14)


# --------------------------------------------------------------------------
# Integration: Surface.reparameterize_theta
# --------------------------------------------------------------------------


def _plasma():
    return PlasmaSurface.from_wout(str(EQUILIBRIA / "wout_d23p4_tm.nc"), ntheta=64, nzeta=64)


@pytest.mark.parametrize("scheme", ["uniform_arclength", "curvature"])
def test_reparameterize_theta_preserves_the_physical_surface(scheme):
    """Area and volume are properties of the shape, not the parameterization.

    The tolerance is set by the *refit*, not by the map: composing a Fourier
    surface with a theta-map generally widens its spectrum (the composition of
    a band-limited function with a non-affine map is not band-limited), so the
    resampled surface needs more modes than the original did. See
    `test_reparameterize_theta_converges_to_the_original_surface`.
    """
    surface = _plasma()
    reparameterized = surface.reparameterize_theta(scheme, ntheta=128, nzeta=64)

    np.testing.assert_allclose(reparameterized.area, surface.area, rtol=1e-2)
    np.testing.assert_allclose(reparameterized.volume, surface.volume, rtol=1e-2)
    assert reparameterized.nfp == surface.nfp
    assert reparameterized.stellarator_symmetric == surface.stellarator_symmetric
    assert reparameterized.standard_toroidal_angle == surface.standard_toroidal_angle


def test_reparameterize_theta_converges_to_the_original_surface():
    """The residual difference above is truncation, so it must shrink with
    resolution -- a fixed offset would mean the map moved the surface."""
    surface = _plasma()
    errors = [
        abs(surface.reparameterize_theta("uniform_arclength", ntheta=n, nzeta=64).area
            - surface.area) / surface.area
        for n in (64, 128, 256)
    ]
    assert errors[0] > errors[1] > errors[2]
    assert errors[-1] < 1e-4


def test_reparameterize_theta_preserves_stellarator_symmetry():
    reparameterized = _plasma().reparameterize_theta("uniform_arclength")
    assert reparameterized.stellarator_symmetric
    assert np.abs(reparameterized.rmns).max() == 0.0
    assert np.abs(reparameterized.zmnc).max() == 0.0


def test_reparameterize_theta_equalizes_arclength_on_a_real_surface():
    """The map's defining property, measured on the *base* surface at the
    mapped angles -- no refit and no derivative of the map involved, so the
    only error is the O(1/N^2) chord-versus-arc one.
    """
    surface = _plasma()
    zeta = (2 * np.pi / surface.nfp) * np.arange(4) / 4

    def chord_spread(n, theta_old):
        pts = surface.evaluate_at(
            np.ravel(theta_old), np.ravel(np.broadcast_to(zeta, theta_old.shape))
        )["r"].reshape((3,) + theta_old.shape)
        dl = np.sqrt(np.sum((np.roll(pts, -1, axis=1) - pts) ** 2, axis=0))
        return ((dl.max(axis=0) - dl.min(axis=0)) / dl.mean(axis=0)).max()

    spreads = []
    for n in (256, 512, 1024):
        tmap = theta_map(
            lambda th, ze: surface._evaluate(th, ze)["r"],
            UniformArclength(), nfp=surface.nfp, ntheta=n, nzeta=4,
        )
        theta_out = 2 * np.pi * np.arange(n) / n
        spreads.append(chord_spread(n, tmap(theta_out, zeta)))

    # In the original angle the arclength varies by more than 100%.
    theta_out = 2 * np.pi * np.arange(1024) / 1024
    naive = chord_spread(1024, np.broadcast_to(theta_out[:, None], (1024, 4)))
    assert naive > 1.0

    assert spreads[-1] < 2e-3
    # Second-order convergence: quadrupling N cuts the spread by ~4x each time.
    assert spreads[0] / spreads[1] > 2.5
    assert spreads[1] / spreads[2] > 2.5


def test_reparameterize_theta_reports_diagnostics():
    reparameterized = _plasma().reparameterize_theta("uniform_arclength", ntheta=64, nzeta=32)
    diagnostics = reparameterized.theta_map.diagnostics
    assert diagnostics["max_residual"] < 1e-5
    assert diagnostics["dl_spread"].shape == (32,)


def test_reparameterize_theta_requires_paired_evaluation():
    """`Surface.evaluate_at` has no representation-independent default, so a
    surface that does not implement it gets a clear error rather than an
    AttributeError deep inside the resampling."""
    from regcoil.surface import Surface

    class GridOnlySurface(Surface):
        nfp, stellarator_symmetric, ntheta, nzeta = 1, True, 16, 8
        standard_toroidal_angle = True

        def _evaluate(self, theta, zetal):
            R = 5.0 + np.cos(theta)[:, None] + 0 * zetal[None, :]
            Z = np.sin(theta)[:, None] + 0 * zetal[None, :]
            phi = zetal[None, :] + 0 * theta[:, None]
            zero = np.zeros_like(R)
            r = np.stack([R * np.cos(phi), R * np.sin(phi), Z], axis=0)
            return {"r": r, "drdtheta": np.stack([zero] * 3), "drdzeta": np.stack([zero] * 3)}

    with pytest.raises(NotImplementedError, match="evaluate_at"):
        GridOnlySurface().reparameterize_theta("uniform_arclength")


# --------------------------------------------------------------------------
# Integration: CoilSurface.from_uniform_offset(theta_reparameterization=...)
# --------------------------------------------------------------------------


@pytest.mark.parametrize("standard_toroidal_angle", [False, True])
def test_from_uniform_offset_reparameterization_preserves_the_surface(standard_toroidal_angle):
    """The reparameterized offset surface is the same physical surface as the
    plain one -- checked by resolving both well enough that the remaining
    difference is the Fourier truncation, not the parameterization."""
    plasma = PlasmaSurface.from_wout(str(EQUILIBRIA / "wout_d23p4_tm.nc"), ntheta=8, nzeta=8)
    kwargs = dict(
        separation=0.5, ntheta=128, nzeta=128, mpol=32, ntor=32,
        standard_toroidal_angle=standard_toroidal_angle,
    )
    plain = CoilSurface.from_uniform_offset(plasma, theta_reparameterization=None, **kwargs)
    mapped = CoilSurface.from_uniform_offset(
        plasma, theta_reparameterization="uniform_arclength", **kwargs
    )

    np.testing.assert_allclose(mapped.area, plain.area, rtol=1e-3)
    np.testing.assert_allclose(mapped.volume, plain.volume, rtol=3e-3)
    assert mapped.standard_toroidal_angle is standard_toroidal_angle
    assert mapped.theta_map is not None
    assert mapped.theta_map.diagnostics["max_residual"] < 1e-4


def test_from_uniform_offset_without_reparameterization_has_no_theta_map():
    plasma = PlasmaSurface.circular_torus(R0=5.0, a=1.0, nfp=3, ntheta=8, nzeta=8)
    coil = CoilSurface.from_uniform_offset(
        plasma, separation=0.3, ntheta=8, nzeta=8, mpol=2, ntor=1,
        theta_reparameterization=None,
    )
    assert coil.theta_map is None


def test_from_uniform_offset_reparameterizes_by_default():
    """`theta_reparameterization` defaults to `"uniform_arclength"`, so the
    returned surface carries the corresponding `theta_map`."""
    plasma = PlasmaSurface.circular_torus(R0=5.0, a=1.0, nfp=3, ntheta=8, nzeta=8)
    coil = CoilSurface.from_uniform_offset(plasma, separation=0.3, ntheta=8, nzeta=8, mpol=2, ntor=1)
    assert coil.theta_map is not None
    assert coil.theta_map.scheme == UniformArclength()


def test_from_uniform_offset_reparameterization_equalizes_arclength():
    """The point of the feature, on the construction it exists for."""
    plasma = PlasmaSurface.from_wout(str(EQUILIBRIA / "wout_d23p4_tm.nc"), ntheta=8, nzeta=8)
    kwargs = dict(separation=0.5, ntheta=64, nzeta=64, mpol=24, ntor=24)
    plain = CoilSurface.from_uniform_offset(plasma, theta_reparameterization=None, **kwargs)
    mapped = CoilSurface.from_uniform_offset(
        plasma, theta_reparameterization="uniform_arclength", **kwargs
    )

    def spread(surf):
        dr = surf.drdtheta[:, :, : surf.nzeta]
        speed = np.sqrt(np.sum(dr * dr, axis=0))
        return ((speed.max(axis=0) - speed.min(axis=0)) / speed.mean(axis=0)).max()

    assert spread(plain) > 0.5
    assert spread(mapped) < 0.05


def test_from_uniform_offset_reparameterization_preserves_stellarator_symmetry():
    plasma = PlasmaSurface.from_wout(str(EQUILIBRIA / "wout_d23p4_tm.nc"), ntheta=8, nzeta=8)
    coil = CoilSurface.from_uniform_offset(
        plasma, separation=0.5, ntheta=32, nzeta=32, mpol=12, ntor=12,
        theta_reparameterization="uniform_arclength",
    )
    assert coil.stellarator_symmetric
    assert np.abs(coil.rmns).max() == 0.0
    assert np.abs(coil.zmnc).max() == 0.0
    assert np.abs(coil.numnc).max() == 0.0


def test_from_uniform_offset_accepts_scheme_instances():
    plasma = PlasmaSurface.circular_torus(R0=5.0, a=1.0, nfp=3, ntheta=8, nzeta=8)
    coil = CoilSurface.from_uniform_offset(
        plasma, separation=0.3, ntheta=16, nzeta=16, mpol=4, ntor=2,
        theta_reparameterization=CurvatureWeighted(exponent=0.5, floor=0.1, refinement=4),
    )
    assert coil.theta_map.scheme == CurvatureWeighted(exponent=0.5, floor=0.1, refinement=4)
