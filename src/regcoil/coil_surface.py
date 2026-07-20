"""`CoilSurface`: the coil winding surface."""

from __future__ import annotations

import logging
from time import perf_counter

import numpy as np

from ._io import read_nescin
from .fourier_surface import FourierSurface

logger = logging.getLogger(__name__)


class CoilSurface(FourierSurface):
    """Coil-side Fourier filtering, plus `from_uniform_offset`
    (the only constructor that calls Fortran, and only when
    `standard_toroidal_angle=True`).
    """

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
          Fourier-transformed directly -- pure Python/numpy, no Fortran call.
          The local normal generally has a toroidal component, so the moved
          point's actual Cartesian toroidal angle `atan2(y, x)` is not `zeta`;
          that offset `atan2(y, x) - zeta` is captured in the surface's `nu`
          (toroidal-angle-shift) Fourier modes (see `FourierSurface`), so the
          moved-along-normal geometry is represented exactly (up to `mpol`/
          `ntor` truncation). The resulting `CoilSurface.standard_toroidal_angle`
          is `False`, and a constant-`zeta` slice of it is *not* a constant
          physical-toroidal-angle cross section.
        - `True`: legacy `geometry_option_coil=2` behavior. Each grid point
          is a Fortran root solve (`regcoil_uniform_offset_surface`, Phase 7)
          for the point on the offset surface whose Cartesian toroidal angle
          equals the coil `zeta` grid value, so `tol` (root-solve tolerance)
          applies only in this mode. `standard_toroidal_angle=True` is set on
          the result.
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
            from . import _core

            logger.info(
                "Starting uniform offset surface kernel for separation=%s, transform grid=%dx%d, modes=%dx%d",
                separation,
                ntheta_transform,
                nzeta_transform,
                mpol,
                ntor,
            )
            kernel_start = perf_counter()
            xm_out, xn_out, rmnc_out, rmns_out, zmnc_out, zmns_out = _core.uniform_offset_surface(  # type: ignore[attr-defined]
                plasma.xm, plasma.xn, plasma.rmnc, plasma.rmns, plasma.zmnc, plasma.zmns,
                lasym, plasma.nfp,
                float(separation), int(mpol), int(ntor),
                int(ntheta_transform), int(nzeta_transform), float(tol),
            )
            logger.info(
                "Finished uniform offset surface kernel in %.3f s",
                perf_counter() - kernel_start,
            )
            numns_out = numnc_out = None
        else:
            theta_grid, zeta_grid, major_r, z_val, nu_val = _offset_along_normal(
                plasma, separation, int(ntheta_transform), int(nzeta_transform)
            )
            (
                xm_out, xn_out, rmnc_out, rmns_out, zmnc_out, zmns_out, numns_out, numnc_out,
            ) = _fourier_transform_offset_surface(
                theta_grid, zeta_grid, major_r, z_val, nu_val, plasma.nfp, int(mpol), int(ntor), lasym
            )

        return cls(
            xm_out, xn_out, rmnc_out, zmns_out, rmns_out, zmnc_out,
            nfp=plasma.nfp, ntheta=ntheta, nzeta=nzeta,
            standard_toroidal_angle=standard_toroidal_angle,
            numns=numns_out, numnc=numnc_out,
        )

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


def _offset_along_normal(plasma, separation, ntheta_transform, nzeta_transform):
    """Sample `plasma` on a `(theta, zeta)` grid spanning one field period,
    move each point `separation` along its local unit normal, and return the
    grid plus the moved points' `(major_R, Z, nu)`, where
    `nu = atan2(y, x) - zeta` is the deviation of the moved point's physical
    toroidal angle from the surface parameter `zeta` (see `FourierSurface`).

    Relies on `plasma`'s nfp-periodicity (the Fourier representation's `xn`
    is a multiple of `nfp`) to guarantee that the moved points are themselves
    exactly periodic under a `2*pi/nfp` rotation, so one field period is
    enough to reconstruct the whole surface by Fourier transform -- the same
    shortcut the Fortran root-solve kernel relies on. `nu` inherits this
    periodicity because both `atan2(y, x)` and `zeta` advance by `2*pi/nfp`
    across one field period, so their difference is periodic.
    """
    theta_grid = 2 * np.pi * np.arange(ntheta_transform) / ntheta_transform
    zeta_grid = (2 * np.pi / plasma.nfp) * np.arange(nzeta_transform) / nzeta_transform

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
    major_r = np.sqrt(r_offset[0] ** 2 + r_offset[1] ** 2)
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
    ordering/convention of the Fortran `regcoil_uniform_offset_surface`
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
    Fortran `regcoil_uniform_offset_surface` kernel (not performance-critical,
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
