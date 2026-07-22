"""Poloidal-angle reparameterization.

A *theta reparameterization* replaces a surface's poloidal angle by a new one,

    theta_old = g(theta_new, zeta),

with `g` strictly increasing in `theta_new` and `g(theta + 2*pi, zeta) =
g(theta, zeta) + 2*pi`. The physical surface is unchanged -- only the labels
of its points move -- but the Fourier spectrum, the current-potential basis
`sin(m*theta - n*zeta)`, and the coil spacing all change with it. The usual
reason to want one is that a uniform-arclength angle needs far fewer poloidal
modes to resolve the same shape.

Because a theta-map slides points *along constant-`zeta` curves*, both schemes
here are defined on the curve `theta -> r(theta, zeta)` at fixed `zeta`, not on
the constant-`phi` cross section. The two coincide exactly when `zeta` is the
standard toroidal angle (`nu = 0`); when it is not, the constant-`phi` cross
section is not a coordinate curve and no theta-map can control it.

Both schemes are the same quadrature. With

    w(theta_old, zeta) = |dr/dtheta| * kappa^exponent,
    theta_new(theta_old) = 2*pi * int_0^theta_old w / int_0^2pi w,

the incremental arclength in the new angle satisfies

    dl/dtheta_new  proportional to  kappa^-exponent,

so `exponent=0` gives uniform arclength and `exponent=1` gives an arclength
inversely proportional to the curvature of the constant-`zeta` curve. No
fixed-point iteration is needed (unlike the legacy Fortran, which iterated only
because it estimated `dl` by finite differences on the coarse grid).

The core here takes a *curve callable*, not a `Surface`, so the same machinery
serves both an existing surface (`Surface.reparameterize_theta`) and the
not-yet-a-surface offset point cloud inside
`CoilSurface.from_uniform_offset(..., theta_reparameterization=...)`, where the
map must be applied *before* the Fourier fit for the fit to happen in the good
angle.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

__all__ = ["UniformArclength", "CurvatureWeighted", "ThetaMap", "theta_map"]


@dataclass(frozen=True)
class UniformArclength:
    """Make `|dr/dtheta|` at fixed `zeta` independent of `theta`.

    Parameters
    ----------
    measure : {"3d", "poloidal"}
        `"3d"` (default) uses the true arclength of the constant-`zeta` curve,
        `sqrt(R'^2 + Z'^2 + R^2 nu'^2)`. `"poloidal"` uses its projection into
        the R-Z plane, `sqrt(R'^2 + Z'^2)` -- what the legacy Fortran
        `geometry_option_coil=4` computed. The two are identical when `nu = 0`.
    refinement : int
        The quadrature grid is `refinement * ntheta` points per `zeta`.
    """

    measure: str = "3d"
    refinement: int = 8


@dataclass(frozen=True)
class CurvatureWeighted:
    """Make the incremental arclength `dl/dtheta` at fixed `zeta` vary as
    `kappa^-exponent`, where `kappa` is the curvature of the constant-`zeta`
    curve -- so points bunch up where the curve bends sharply.

    Parameters
    ----------
    exponent : float
        `1.0` (default) gives `dl/dtheta` exactly inversely proportional to
        `kappa`; `0.0` degenerates to `UniformArclength`. Fractional values
        blend the two.
    floor : float
        `kappa` is replaced by `kappa + floor * <kappa>` (with `<kappa>` the
        arclength-weighted mean). Without this the scheme is singular: a
        straight segment has `kappa = 0` and would be allotted an unbounded
        share of the arclength.
    measure, refinement
        As in `UniformArclength`.
    """

    exponent: float = 1.0
    floor: float = 0.05
    measure: str = "3d"
    refinement: int = 8


_NAMED_SCHEMES = {
    "uniform_arclength": UniformArclength,
    "curvature": CurvatureWeighted,
}


def _as_scheme(scheme):
    """Accept a scheme dataclass, or the shorthand strings
    `"uniform_arclength"` / `"curvature"`."""
    if isinstance(scheme, (UniformArclength, CurvatureWeighted)):
        return scheme
    if isinstance(scheme, str):
        try:
            return _NAMED_SCHEMES[scheme]()
        except KeyError:
            raise ValueError(
                f"Unknown theta reparameterization {scheme!r}; expected one of "
                f"{sorted(_NAMED_SCHEMES)}, or a UniformArclength/CurvatureWeighted "
                f"instance."
            ) from None
    raise TypeError(
        f"theta reparameterization must be a UniformArclength or CurvatureWeighted "
        f"instance, or one of the strings {sorted(_NAMED_SCHEMES)}; got "
        f"{type(scheme).__name__}."
    )


class ThetaMap:
    """The map `theta_old = g(theta_new, zeta)`, stored as Fourier modes of its
    periodic part `g - theta_new` (which is genuinely periodic in both angles,
    unlike `g` itself, whose secular term is the identity).

    The sine/cosine convention matches `FourierSurface`: the angle is
    `m*theta - n*zeta` and `xn` already includes the factor of `nfp`. Under
    stellarator symmetry `g - theta` is odd, so `gmnc = 0` -- the same parity
    as `Z`'s `zmnc`.
    """

    def __init__(self, xm, xn, gmns, gmnc, *, nfp, scheme=None, diagnostics=None):
        self.xm = np.asarray(xm, dtype=np.int64)
        self.xn = np.asarray(xn, dtype=np.int64)
        self.gmns = np.asarray(gmns, dtype=np.float64)
        self.gmnc = np.asarray(gmnc, dtype=np.float64)
        self.nfp = int(nfp)
        #: The scheme this map was built from, kept for provenance.
        self.scheme = scheme
        #: Quality of the achieved reparameterization; see `theta_map`.
        self.diagnostics = {} if diagnostics is None else diagnostics

    def __call__(self, theta, zeta):
        """Evaluate `theta_old` on the tensor grid `(theta, zeta)`.

        Returns a `(len(theta), len(zeta))` array.
        """
        theta = np.atleast_1d(np.asarray(theta, dtype=float))
        zeta = np.atleast_1d(np.asarray(zeta, dtype=float))
        angle = (
            self.xm[:, None, None] * theta[None, :, None]
            - self.xn[:, None, None] * zeta[None, None, :]
        )
        periodic = np.tensordot(self.gmns, np.sin(angle), axes=(0, 0)) + np.tensordot(
            self.gmnc, np.cos(angle), axes=(0, 0)
        )
        return theta[:, None] + periodic


def _fft_derivative(field, order, axis=0):
    """`d^order field / dtheta^order` for a field sampled on a uniform periodic
    theta grid (`axis`), by spectral differentiation.

    Positions are all we ever require of a curve this way:
    the alternative, analytic derivatives, would need the plasma surface's
    third derivatives on the `from_uniform_offset` path, and would miss the
    implicit-function correction on the `standard_toroidal_angle=True` path,
    where the curve is defined by a per-point root solve.
    """
    n = field.shape[axis]
    k = np.fft.fftfreq(n, d=1.0 / n)
    factor = (1j * k) ** order
    if order % 2 == 1 and n % 2 == 0:
        # The Nyquist mode carries no meaningful odd derivative for a real
        # signal (numpy assigns it k = -n/2); zeroing it keeps the result real.
        factor[n // 2] = 0.0
    shape = [1] * field.ndim
    shape[axis] = n
    transformed = np.fft.fft(field, axis=axis) * factor.reshape(shape)
    return np.real(np.fft.ifft(transformed, axis=axis))


def _cumulative_integral(w):
    """`F_j = int_0^theta_j w dtheta` for `w` sampled on a uniform periodic
    theta grid (axis 0), by spectral integration.

    Writing `w = c_0 + sum_{k != 0} c_k e^{i k theta}`, the antiderivative is
    `c_0*theta + sum_{k != 0} c_k (e^{i k theta} - 1) / (i k)`. This is
    spectrally accurate, unlike the O(h^2) cumulative trapezoid rule the legacy
    code effectively used.
    """
    n = w.shape[0]
    theta = 2 * np.pi * np.arange(n) / n
    k = np.fft.fftfreq(n, d=1.0 / n)
    c = np.fft.fft(w, axis=0) / n

    shape = [1] * w.ndim
    shape[0] = n
    nonzero = k != 0
    inverse = np.zeros(n, dtype=complex)
    inverse[nonzero] = 1.0 / (1j * k[nonzero])
    d = c * inverse.reshape(shape)

    # `n * ifft(d)` evaluates sum_k d_k e^{i k theta_j}; subtracting its value
    # at theta = 0 (which is sum_k d_k) imposes F(0) = 0.
    oscillatory = n * np.fft.ifft(d, axis=0) - np.sum(d, axis=0)
    result = np.real(oscillatory) + np.real(c[0]) * theta.reshape(shape)
    # That cancellation is only exact in exact arithmetic; at large `n` it
    # leaves F[0] at roundoff rather than zero. Subtracting it makes the
    # anchor F(0) = 0 exact, which the periodic spline that inverts this map
    # requires to machine precision (and which is what preserves stellarator
    # symmetry).
    return result - result[0]


def _weight(positions, scheme):
    """The quadrature weight `w = |dr/dtheta| * kappa^exponent` on the fine
    theta grid, plus the incremental arclength `|dr/dtheta|` itself.

    `positions` is `(3, ntheta_fine, nzeta)`. For `measure="poloidal"` the
    curve is first projected to `(R, Z)`, so both speed and curvature are those
    of the R-Z cross-section curve rather than of the 3D constant-`zeta` curve.
    """
    if scheme.measure == "poloidal":
        R = np.hypot(positions[0], positions[1])
        curve = np.stack([R, positions[2]], axis=0)
    elif scheme.measure == "3d":
        curve = positions
    else:
        raise ValueError(f"measure must be '3d' or 'poloidal', got {scheme.measure!r}")

    d1 = _fft_derivative(curve, 1, axis=1)
    speed = np.sqrt(np.sum(d1 * d1, axis=0))

    exponent = float(getattr(scheme, "exponent", 0.0))
    if exponent == 0.0:
        return speed, speed

    d2 = _fft_derivative(curve, 2, axis=1)
    # |r' x r''| / |r'|^3, written so it covers the 2-component (poloidal)
    # case as well as the 3-component one.
    if curve.shape[0] == 2:
        cross = np.abs(d1[0] * d2[1] - d1[1] * d2[0])
    else:
        cross = np.sqrt(
            (d1[1] * d2[2] - d1[2] * d2[1]) ** 2
            + (d1[2] * d2[0] - d1[0] * d2[2]) ** 2
            + (d1[0] * d2[1] - d1[1] * d2[0]) ** 2
        )
    kappa = cross / np.maximum(speed, np.finfo(float).tiny) ** 3

    # Arclength-weighted mean curvature, per zeta, sets the floor's scale.
    mean_kappa = np.sum(kappa * speed, axis=0) / np.sum(speed, axis=0)
    kappa = kappa + float(scheme.floor) * mean_kappa[None, :]
    if np.any(kappa <= 0):
        raise ValueError(
            "Curvature-weighted reparameterization requires a positive weight "
            "everywhere, but the floored curvature reached zero. Increase "
            "`CurvatureWeighted.floor`."
        )
    return speed * kappa**exponent, speed


def theta_map(curve, scheme, *, nfp, ntheta, nzeta, mpol=None, ntor=None,
              stellarator_symmetric=True):
    """Build the `ThetaMap` that reparameterizes `curve`'s poloidal angle
    according to `scheme`.

    Parameters
    ----------
    curve : callable
        `curve(theta, zeta) -> (3, len(theta), len(zeta))` Cartesian positions
        on a tensor grid. This is the only thing the machinery requires of a
        geometry, which is what lets it serve both an existing `Surface` and a
        point cloud that is not a surface yet.
    scheme : UniformArclength | CurvatureWeighted | str
    nfp, ntheta, nzeta : int
        Field periods, and the output grid the map is tabulated on. `zeta`
        spans one field period.
    mpol, ntor : int, optional
        Mode numbers for the stored map. Default to the output grid's Nyquist
        limit, at which the transform round-trips exactly, so the default adds
        no truncation error. Lower values deliberately smooth the map -- mainly
        useful in `zeta`, since each `zeta` is solved independently and nothing
        otherwise couples them.
    stellarator_symmetric : bool
        When True the cosine modes of `g - theta` are zeroed, enforcing the
        parity that the anchor `g(0, zeta) = 0` already produces analytically.

    Returns
    -------
    ThetaMap, whose `diagnostics` dict reports the achieved quality; see
    `_diagnostics` for `dl_spread` (what the legacy
    `constant_arclength_tolerance` measured) and `residual` (the scheme-neutral
    correctness check).
    """
    from scipy.interpolate import CubicSpline

    scheme = _as_scheme(scheme)
    refinement = int(scheme.refinement)
    if refinement < 1:
        raise ValueError(f"refinement must be >= 1, got {refinement}")

    ntheta_fine = refinement * int(ntheta)
    theta_fine = 2 * np.pi * np.arange(ntheta_fine) / ntheta_fine
    zeta_grid = (2 * np.pi / nfp) * np.arange(nzeta) / nzeta
    theta_out = 2 * np.pi * np.arange(ntheta) / ntheta

    positions = np.asarray(curve(theta_fine, zeta_grid))
    if positions.shape != (3, ntheta_fine, nzeta):
        raise ValueError(
            f"curve(theta, zeta) must return shape (3, {ntheta_fine}, {nzeta}), "
            f"got {positions.shape}."
        )

    w, speed = _weight(positions, scheme)
    if np.any(w <= 0):
        raise ValueError(
            "Theta reparameterization requires a positive weight everywhere, but "
            "|dr/dtheta| vanishes somewhere on the surface (a degenerate "
            "constant-zeta curve)."
        )

    # theta_new as a function of theta_old, anchored at theta_new(0) = 0 so
    # that stellarator symmetry survives the reparameterization.
    # The full-turn integral is exactly `mean(w) * 2*pi`, so normalizing by the
    # mean is both exact and free of any endpoint special-casing.
    integral = _cumulative_integral(w)
    theta_new = integral / w.mean(axis=0)[None, :]

    # Invert. The interpolated quantity is theta_old - theta_new, which is
    # periodic (both advance by 2*pi over a poloidal turn) -- the same trick
    # the legacy Fortran used with its periodic spline.
    theta_old = np.empty((int(ntheta), int(nzeta)))
    for j in range(int(nzeta)):
        knots = np.concatenate([theta_new[:, j], [2 * np.pi]])
        values = np.concatenate([theta_fine - theta_new[:, j], [0.0]])
        spline = CubicSpline(knots, values, bc_type="periodic")
        theta_old[:, j] = spline(theta_out) + theta_out

    diagnostics = _diagnostics(theta_new, theta_fine, speed, w, int(ntheta))

    mpol = int(ntheta) // 2 if mpol is None else int(mpol)
    ntor = int(nzeta) // 2 if ntor is None else int(ntor)
    xm, xn, gmns, gmnc = _fit_map(
        theta_old - theta_out[:, None], theta_out, zeta_grid, nfp, mpol, ntor,
        stellarator_symmetric,
    )
    return ThetaMap(xm, xn, gmns, gmnc, nfp=nfp, scheme=scheme, diagnostics=diagnostics)


def _diagnostics(theta_new, theta_fine, speed, w, ntheta):
    """How well the map achieved what the scheme asked for.

    `dl_spread` is the relative spread of the incremental arclength
    `dl/dtheta_new` -- the quantity the legacy `constant_arclength_tolerance`
    thresholded, reported here rather than used to control a loop.
    It is near zero for `UniformArclength` but *deliberately*
    large for `CurvatureWeighted`, which asks for a non-uniform `dl`.

    `residual` is the scheme-neutral correctness check: by construction
    `w / (dtheta_new/dtheta_old)` equals the constant `int w / (2*pi)` for
    *any* scheme, so its relative spread is the map's error. For
    `UniformArclength` the two coincide.
    """
    dtheta_new = _fft_derivative(theta_new - theta_fine[:, None], 1) + 1.0
    dl = speed / dtheta_new
    invariant = w / dtheta_new
    dl_spread = (dl.max(axis=0) - dl.min(axis=0)) / dl.mean(axis=0)
    residual = (invariant.max(axis=0) - invariant.min(axis=0)) / invariant.mean(axis=0)
    return {
        "dl_spread": dl_spread,
        "max_dl_spread": float(dl_spread.max()),
        "residual": residual,
        "max_residual": float(residual.max()),
        "ntheta": ntheta,
    }


def _fit_map(periodic_part, theta_grid, zeta_grid, nfp, mpol, ntor, stellarator_symmetric):
    """DFT the periodic part of the map onto `mpol`/`ntor` modes, with the same
    convention and Nyquist weighting as the offset-surface transform."""
    from .coil_surface import _uniform_offset_modes

    ntheta = theta_grid.shape[0]
    nzeta = zeta_grid.shape[0]
    xm, xn = _uniform_offset_modes(nfp, mpol, ntor)

    angle = (
        xm[:, None, None] * theta_grid[None, :, None]
        - xn[:, None, None] * zeta_grid[None, None, :]
    )
    weight = np.full(xm.shape[0], 2.0 / (ntheta * nzeta))
    if ntheta % 2 == 0:
        weight[xm == ntheta // 2] /= 2
    if nzeta % 2 == 0:
        weight[np.abs(xn) == nfp * (nzeta // 2)] /= 2

    gmnc = np.tensordot(np.cos(angle), periodic_part, axes=([1, 2], [0, 1])) * weight
    gmns = np.tensordot(np.sin(angle), periodic_part, axes=([1, 2], [0, 1])) * weight
    gmnc[0] = periodic_part.mean()
    gmns[0] = 0.0
    if stellarator_symmetric:
        gmnc[:] = 0.0
    return xm, xn, gmns, gmnc
