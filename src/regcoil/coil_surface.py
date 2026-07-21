"""`CoilSurface`: the coil winding surface."""

from __future__ import annotations

import logging
from time import perf_counter

import numpy as np

from ._io import read_nescin
from .fourier_surface import FourierSurface

logger = logging.getLogger(__name__)


class CoilSurface(FourierSurface):
    """Coil-side Fourier filtering, plus `from_uniform_offset`."""

    @classmethod
    def from_nescin(cls, filename, nfp, ntheta=64, nzeta=64, mpol_filter=None, ntor_filter=None):
        """Read a nescin-format coil surface. `nfp` must match the plasma
        surface's `nfp` (a nescin file does not encode it).
        """
        data = read_nescin(filename, nfp)
        surf = cls(
            data["xm"], data["xn"], data["rmnc"], data["zmns"], data["rmns"], data["zmnc"],
            nfp=nfp, ntheta=ntheta, nzeta=nzeta,
        )
        if mpol_filter is not None or ntor_filter is not None:
            surf.filter_modes(mpol_filter, ntor_filter)
        return surf

    @classmethod
    def from_uniform_offset(
        cls, plasma, separation, ntheta=64, nzeta=64, mpol=24, ntor=24,
        standard_toroidal_angle=False,
        ntheta_transform=None, nzeta_transform=None, tol=1e-10,
        theta_reparameterization="uniform_arclength",
    ):
        """A surface offset uniformly outward from `plasma` by `separation`
        meters along its normal.

        `ntheta_transform`/`nzeta_transform` (default `ntheta`/`nzeta`) set
        the grid the offset surface is sampled on before being
        Fourier-transformed to `mpol`/`ntor` modes.

        `standard_toroidal_angle` selects the construction:

        - `False` (the default): each `plasma` grid point at `(theta, zeta)`
          is moved `separation` along the local unit normal, and the moved
          points are labeled by that same `(theta, zeta)` and
          Fourier-transformed directly -- no root solve is needed.
          The local normal generally has a toroidal component, so the moved
          point's actual Cartesian toroidal angle `atan2(y, x)` is not `zeta`;
          that offset `atan2(y, x) - zeta` is captured in the surface's `nu`
          (toroidal-angle-shift) Fourier modes (see `FourierSurface`), so the
          moved-along-normal geometry is represented exactly (up to `mpol`/
          `ntor` truncation). The resulting `CoilSurface.standard_toroidal_angle`
          is `False`, and a constant-`zeta` slice of it is *not* a constant
          physical-toroidal-angle cross section.
        - `True`: legacy `geometry_option_coil=2` behavior. Each grid point is
          a root solve (`_solve_zeta_for_phi`) for the point on the offset
          surface whose Cartesian toroidal angle equals the coil `zeta` grid
          value, so `tol` (root-solve tolerance) applies only in this mode.
          `standard_toroidal_angle=True` is set on the result. Raises if the
          offset surface is self-intersecting, which leaves some target angles
          unreachable.

        `theta_reparameterization` (default `"uniform_arclength"`)
        reparameterizes the coil surface's poloidal angle without changing its
        shape -- a `UniformArclength` or `CurvatureWeighted` instance, or the
        shorthand string `"uniform_arclength"` / `"curvature"` (see
        `regcoil.reparameterize`); pass `None` to keep the plasma's poloidal
        angle. The legacy `geometry_option_coil=2` is
        `standard_toroidal_angle=True, theta_reparameterization=None`, and
        `geometry_option_coil=4` is
        `standard_toroidal_angle=True, theta_reparameterization="uniform_arclength"`.
        The map is built from, and applied to, the offset points *before* the
        Fourier transform, so the fit happens once and in the good angle; the
        resulting `ThetaMap` is kept on the returned surface as `theta_map`.
        """
        lasym = not plasma.stellarator_symmetric
        ntheta_transform = ntheta if ntheta_transform is None else ntheta_transform
        nzeta_transform = nzeta if nzeta_transform is None else nzeta_transform

        if standard_toroidal_angle:
            if not plasma.standard_toroidal_angle:
                raise ValueError(
                    "Cannot construct a standard-toroidal-angle uniform-offset surface from a "
                    "non-standard-toroidal-angle plasma surface. Use "
                    "`CoilSurface.from_uniform_offset(..., standard_toroidal_angle=False)` "
                    "instead."
                )
            if not isinstance(plasma, FourierSurface):
                raise ValueError(
                    "Cannot construct a standard-toroidal-angle uniform-offset surface from a "
                    "non-Fourier plasma surface. Use "
                    "`CoilSurface.from_uniform_offset(..., standard_toroidal_angle=False)` "
                    "instead."
                )
            logger.info(
                "Starting uniform offset surface root solve for separation=%s, transform grid=%dx%d, modes=%dx%d",
                separation,
                ntheta_transform,
                nzeta_transform,
                mpol,
                ntor,
            )
            kernel_start = perf_counter()

            def curve(theta, zeta):
                theta_2d, zeta_2d = np.meshgrid(theta, zeta, indexing="ij")
                return _offset_curve_at_standard_toroidal_angle(
                    plasma, separation, theta_2d, zeta_2d, float(tol)
                )

            tmap = _build_theta_map(
                theta_reparameterization, curve, plasma,
                int(ntheta_transform), int(nzeta_transform),
            )
            theta_grid, zeta_grid, major_r, z_val = _offset_at_standard_toroidal_angle(
                plasma, separation, int(ntheta_transform), int(nzeta_transform), float(tol),
                theta=_map_theta(tmap, int(ntheta_transform), int(nzeta_transform), plasma.nfp),
            )
            logger.info(
                "Finished uniform offset surface root solve in %.3f s",
                perf_counter() - kernel_start,
            )
            nu_val = np.zeros_like(major_r)
        else:

            def curve(theta, zeta):
                theta_2d, zeta_2d = np.meshgrid(theta, zeta, indexing="ij")
                return _offset_points(plasma, separation, theta_2d, zeta_2d)

            tmap = _build_theta_map(
                theta_reparameterization, curve, plasma,
                int(ntheta_transform), int(nzeta_transform),
            )
            theta_grid, zeta_grid, major_r, z_val, nu_val = _offset_along_normal(
                plasma, separation, int(ntheta_transform), int(nzeta_transform),
                theta=_map_theta(tmap, int(ntheta_transform), int(nzeta_transform), plasma.nfp),
            )

        (
            xm_out, xn_out, rmnc_out, rmns_out, zmnc_out, zmns_out, numns_out, numnc_out,
        ) = _fourier_transform_offset_surface(
            theta_grid, zeta_grid, major_r, z_val, nu_val, plasma.nfp, int(mpol), int(ntor), lasym
        )
        if standard_toroidal_angle:
            # phi == zeta by construction, so nu is identically zero (up to the
            # root-solve tolerance); do not carry meaningless modes.
            numns_out = numnc_out = None

        surface = cls(
            xm_out, xn_out, rmnc_out, zmns_out, rmns_out, zmnc_out,
            nfp=plasma.nfp, ntheta=ntheta, nzeta=nzeta,
            standard_toroidal_angle=standard_toroidal_angle,
            numns=numns_out, numnc=numnc_out,
        )
        surface.theta_map = tmap
        return surface

    def filter_modes(self, mpol_filter=None, ntor_filter=None):
        """Zero out modes with `|m|` > mpol_filter or `|n|` > ntor_filter*nfp, in place.

        Filters the toroidal-angle-shift modes (`numns`/`numnc`) alongside
        `R`/`Z`, so a non-standard-toroidal-angle surface stays consistent.
        """
        mpol_limit = np.inf if mpol_filter is None else mpol_filter
        ntor_limit = np.inf if ntor_filter is None else ntor_filter * self.nfp
        mask = (np.abs(self.xm) > mpol_limit) | (np.abs(self.xn) > ntor_limit)
        self.rmnc[mask] = 0
        self.rmns[mask] = 0
        self.zmnc[mask] = 0
        self.zmns[mask] = 0
        self.numns[mask] = 0
        self.numnc[mask] = 0

    def save(self, path):
        """Save this coil surface to `path` (NetCDF-4 via `h5netcdf`; ADR-028)."""
        from . import _serialize
        _serialize.save(path, coil=self)

    @classmethod
    def load(cls, path):
        """Load a `CoilSurface` from `path`."""
        from . import _serialize
        return _serialize.load_coil(path)


def _build_theta_map(scheme, curve, plasma, ntheta_transform, nzeta_transform):
    """`None` if no reparameterization was requested, else the `ThetaMap` for
    `curve`. Split out only so the two `from_uniform_offset` branches share it."""
    if scheme is None:
        return None
    from .reparameterize import theta_map

    return theta_map(
        curve, scheme, nfp=plasma.nfp,
        ntheta=ntheta_transform, nzeta=nzeta_transform,
        stellarator_symmetric=plasma.stellarator_symmetric,
    )


def _map_theta(tmap, ntheta_transform, nzeta_transform, nfp):
    """The `(ntheta_transform, nzeta_transform)` array of plasma poloidal
    angles to sample at, or `None` for the plain uniform grid."""
    if tmap is None:
        return None
    theta_grid = 2 * np.pi * np.arange(ntheta_transform) / ntheta_transform
    zeta_grid = (2 * np.pi / nfp) * np.arange(nzeta_transform) / nzeta_transform
    return tmap(theta_grid, zeta_grid)


def _offset_points(plasma, separation, theta, zeta):
    """Move `plasma` a distance `separation` along its local unit normal, at
    *paired* `(theta, zeta)` points of any matching shape (unlike
    `_offset_along_normal`, which works on a tensor-product grid).

    Returns a `(3,) + theta.shape` array of Cartesian positions.
    """
    shape = np.shape(theta)
    evaluated = plasma.evaluate_at(np.ravel(theta), np.ravel(zeta))
    r = evaluated["r"]
    drdtheta = evaluated["drdtheta"]
    drdzeta = evaluated["drdzeta"]

    # N = dr/dzeta x dr/dtheta, matching Surface.normal's sign convention.
    normal = np.empty_like(r)
    normal[0] = drdzeta[1] * drdtheta[2] - drdtheta[1] * drdzeta[2]
    normal[1] = drdzeta[2] * drdtheta[0] - drdtheta[2] * drdzeta[0]
    normal[2] = drdzeta[0] * drdtheta[1] - drdtheta[0] * drdzeta[1]
    normal /= np.sqrt(np.sum(normal * normal, axis=0))

    return (r + separation * normal).reshape((3,) + shape)


def _solve_zeta_for_phi(plasma, separation, theta, phi_target, tol, max_iterations=100):
    """For each `(theta, phi_target)`, find the surface parameter `zeta` such
    that the offset point at `(theta, zeta)` has Cartesian toroidal angle
    `atan2(y, x) = phi_target` -- the legacy "cosm" construction that makes a
    constant-`zeta` slice a constant-physical-angle cross section (ADR-031
    decision 2; replaces the Fortran `regcoil_fzero` root solve).

    Vectorized bisection-safeguarded secant on the bracket
    `[phi_target - 1, phi_target + 1]` (the same bracket the Fortran kernel
    used). A plain fixed-point iteration `zeta <- zeta - residual` is *not*
    robust enough: it needs `0 < dphi/dzeta < 2`, and while a large-aspect-ratio
    case like `d23p4_tm` at `separation=0.5` stays within `[0.80, 1.24]`, a
    compact device such as `li383_1.4m` at `separation/a ~ 1.5` reaches 2.7 and
    the iteration diverges. Bracketing also gives us the Fortran's `info=3`
    self-intersecting-offset-surface guard for free, as a missing sign change.
    """
    phi_target = np.asarray(phi_target, dtype=float)

    def residual(zeta):
        r = _offset_points(plasma, separation, theta, zeta)
        return np.mod(np.arctan2(r[1], r[0]) - phi_target + np.pi, 2 * np.pi) - np.pi

    a = phi_target - 1.0
    b = phi_target + 1.0
    fa = residual(a)
    fb = residual(b)

    bad = (fa * fb) > 0
    if np.any(bad):
        raise ValueError(
            f"Cannot construct a standard-toroidal-angle uniform-offset surface with "
            f"separation={separation}: no solution for the toroidal angle at "
            f"{bad.sum()} of {bad.size} grid points (the residual does not change sign "
            f"within +/-1 radian). The offset surface is probably self-intersecting -- "
            f"try a smaller separation, or "
            f"`CoilSurface.from_uniform_offset(..., standard_toroidal_angle=False)`, "
            f"which needs no root solve."
        )

    # Illinois-modified regula falsi. Plain false position stagnates here --
    # `phi(zeta)` is convex enough over a 1-radian bracket that one endpoint
    # sticks and convergence degrades to linear (~60 residual evaluations
    # measured). Halving the retained endpoint's function value restores
    # superlinear convergence (~10 evaluations) while never leaving the
    # bracket, so the self-intersecting guard above still holds throughout.
    c, fc = b, fb
    for _ in range(max_iterations):
        # Converge on the residual rather than the bracket width: `dphi/dzeta`
        # is bounded well away from zero (~0.8 to ~2.7 across the cases we
        # test), so |residual| < tol bounds the error in `zeta` by about the
        # same tol, and it is reached in far fewer steps than the ~34
        # bisections a bracket-width criterion would force.
        if np.all(np.abs(fc) < tol):
            break
        with np.errstate(divide="ignore", invalid="ignore"):
            c = b - fb * (b - a) / (fb - fa)
        # Degenerate denominator (both endpoints already at the root): bisect.
        c = np.where(np.isfinite(c), c, 0.5 * (a + b))
        fc = residual(c)

        opposite = (fc * fb) < 0
        a = np.where(opposite, b, a)
        fa = np.where(opposite, fb, 0.5 * fa)
        b, fb = c, fc
    else:
        raise RuntimeError(
            f"Toroidal-angle root solve did not converge to tol={tol} in "
            f"{max_iterations} iterations (max |residual| "
            f"{np.abs(fc).max():.3e})."
        )

    return c


def _offset_at_standard_toroidal_angle(
    plasma, separation, ntheta_transform, nzeta_transform, tol, theta=None
):
    """Sample the uniform-offset surface on a `(theta, zeta)` grid in which
    `zeta` *is* the standard toroidal angle, by solving
    `atan2(y, x) = zeta` at each grid point (`_solve_zeta_for_phi`).

    `theta` (default: the uniform grid) is the `(ntheta_transform,
    nzeta_transform)` array of *plasma* poloidal angles to sample at -- how a
    theta reparameterization enters, since it makes the sampled theta depend
    on zeta.

    Returns `(theta_grid, zeta_grid, major_R, Z)`; there is no `nu` because
    `phi = zeta` holds by construction.
    """
    theta_grid = 2 * np.pi * np.arange(ntheta_transform) / ntheta_transform
    zeta_grid = (2 * np.pi / plasma.nfp) * np.arange(nzeta_transform) / nzeta_transform
    theta_2d = (
        np.broadcast_to(theta_grid[:, None], (ntheta_transform, nzeta_transform))
        if theta is None
        else np.asarray(theta, dtype=float)
    )
    zeta_2d = np.broadcast_to(zeta_grid, theta_2d.shape)

    r_offset = _offset_curve_at_standard_toroidal_angle(
        plasma, separation, theta_2d, zeta_2d, tol
    )
    major_r = np.hypot(r_offset[0], r_offset[1])
    return theta_grid, zeta_grid, major_r, r_offset[2]


def _offset_curve_at_standard_toroidal_angle(plasma, separation, theta, phi_target, tol):
    """The offset point at poloidal angle `theta` whose *physical* toroidal
    angle is `phi_target`: solve for the surface parameter `zeta`, then place
    the point. Paired arrays of any matching shape.

    Kept separate from `_offset_at_standard_toroidal_angle` so that the theta
    reparameterization can use it as its `curve` callable on a refined theta
    grid, before any Fourier fit happens.
    """
    zeta_solution = _solve_zeta_for_phi(plasma, separation, theta, phi_target, tol)
    return _offset_points(plasma, separation, theta, zeta_solution)


def _offset_along_normal(plasma, separation, ntheta_transform, nzeta_transform, theta=None):
    """Sample `plasma` on a `(theta, zeta)` grid spanning one field period,
    move each point `separation` along its local unit normal, and return the
    grid plus the moved points' `(major_R, Z, nu)`, where
    `nu = atan2(y, x) - zeta` is the deviation of the moved point's physical
    toroidal angle from the surface parameter `zeta` (see `FourierSurface`).

    `theta` (default: the uniform grid) is the `(ntheta_transform,
    nzeta_transform)` array of plasma poloidal angles to sample at, which is
    how a theta reparameterization enters. Supplying it requires `plasma` to
    implement `evaluate_at`, since the sample points no longer form a grid.

    Relies on `plasma`'s nfp-periodicity (the Fourier representation's `xn`
    is a multiple of `nfp`) to guarantee that the moved points are themselves
    exactly periodic under a `2*pi/nfp` rotation, so one field period is
    enough to reconstruct the whole surface by Fourier transform -- the same
    shortcut the root-solve construction relies on. `nu` inherits this
    periodicity because both `atan2(y, x)` and `zeta` advance by `2*pi/nfp`
    across one field period, so their difference is periodic.
    """
    theta_grid = 2 * np.pi * np.arange(ntheta_transform) / ntheta_transform
    zeta_grid = (2 * np.pi / plasma.nfp) * np.arange(nzeta_transform) / nzeta_transform

    if theta is None:
        # The sample points form a tensor grid, so use the (faster) grid
        # evaluation, which every Surface supports.
        evaluated = plasma._evaluate(theta_grid, zeta_grid)
        r = evaluated["r"]
        drdtheta = evaluated["drdtheta"]
        drdzeta = evaluated["drdzeta"]

        # N = dr/dzeta x dr/dtheta, matching Surface.normal's sign convention.
        normal = np.empty_like(r)
        normal[0] = drdzeta[1] * drdtheta[2] - drdtheta[1] * drdzeta[2]
        normal[1] = drdzeta[2] * drdtheta[0] - drdtheta[2] * drdzeta[0]
        normal[2] = drdzeta[0] * drdtheta[1] - drdtheta[0] * drdzeta[1]
        normal /= np.sqrt(np.sum(normal * normal, axis=0))

        r_offset = r + separation * normal
    else:
        theta_2d = np.asarray(theta, dtype=float)
        zeta_2d = np.broadcast_to(zeta_grid, theta_2d.shape)
        r_offset = _offset_points(plasma, separation, theta_2d, zeta_2d)

    major_r = np.hypot(r_offset[0], r_offset[1])
    z_val = r_offset[2]

    # Toroidal-angle shift of the moved point relative to zeta. Wrap into
    # (-pi, pi]; the true shift is small (~ separation/R << pi) for any
    # sensible offset, so this picks the correct branch even where atan2
    # crosses its +/-pi cut.
    nu_val = np.arctan2(r_offset[1], r_offset[0]) - zeta_grid[None, :]
    nu_val = np.mod(nu_val + np.pi, 2 * np.pi) - np.pi
    return theta_grid, zeta_grid, major_r, z_val, nu_val


def _uniform_offset_modes(nfp, mpol, ntor):
    """Mode list (`xm`, `xn`) for a uniform-offset surface, matching the
    ordering/convention of the legacy Fortran `regcoil_uniform_offset_surface`
    kernel: `xm >= 0`, `xn` any sign, and `xn` scaled by `nfp` (VMEC
    convention).
    """
    xm = [0]
    xn = [0]
    for jn in range(1, ntor + 1):
        xm.append(0)
        xn.append(jn)
    for jm in range(1, mpol + 1):
        for jn in range(-ntor, ntor + 1):
            xm.append(jm)
            xn.append(jn)
    return np.array(xm, dtype=np.int64), np.array(xn, dtype=np.int64) * nfp


def _fourier_transform_offset_surface(theta_grid, zeta_grid, major_r, z_val, nu_val, nfp, mpol, ntor, lasym):
    """DFT `major_r`/`z_val`/`nu_val` (each `(ntheta_transform,
    nzeta_transform)`, one field period, sampled at `theta_grid`/`zeta_grid`)
    onto `mpol`/`ntor` Fourier modes, with the same normalization as the
    legacy Fortran `regcoil_uniform_offset_surface` kernel (not performance-critical,
    so plain trigonometric sums rather than an FFT, matching the legacy code).

    `nu_val` (the toroidal-angle shift) transforms with the same sin/cos
    convention as `Z`: `numns` is the sine (stellarator-even) part, `numnc`
    the cosine part (zeroed for a stellarator-symmetric plasma).
    """
    ntheta_transform = theta_grid.shape[0]
    nzeta_transform = zeta_grid.shape[0]
    xm, xn = _uniform_offset_modes(nfp, mpol, ntor)
    mnmax = xm.shape[0]

    angle = xm[:, None, None] * theta_grid[None, :, None] - xn[:, None, None] * zeta_grid[None, None, :]
    cosangle = np.cos(angle)
    sinangle = np.sin(angle)

    weight = np.full(mnmax, 2.0 / (ntheta_transform * nzeta_transform))
    # Halve the weight of the Nyquist mode so inverse-transform(transform(.)) is the identity.
    if ntheta_transform % 2 == 0:
        weight[xm == ntheta_transform // 2] /= 2
    if nzeta_transform % 2 == 0:
        weight[np.abs(xn) == nfp * (nzeta_transform // 2)] /= 2

    def _dft(field):
        cmn = np.zeros(mnmax)
        smn = np.zeros(mnmax)
        cmn[1:] = np.tensordot(cosangle[1:], field, axes=([1, 2], [0, 1])) * weight[1:]
        smn[1:] = np.tensordot(sinangle[1:], field, axes=([1, 2], [0, 1])) * weight[1:]
        cmn[0] = field.mean()  # m=n=0: DC term is the plain average, not 2x.
        return cmn, smn

    rmnc, rmns = _dft(major_r)
    zmnc, zmns = _dft(z_val)
    numnc, numns = _dft(nu_val)

    if not lasym:
        # Stellarator symmetry: R even (rmns=0), Z odd (zmnc=0), nu odd (numnc=0).
        rmns[:] = 0
        zmnc[:] = 0
        numnc[:] = 0

    return xm, xn, rmnc, rmns, zmnc, zmns, numns, numnc
