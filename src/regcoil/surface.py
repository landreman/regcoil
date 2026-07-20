"""The `Surface` abstract base class.

The contract is evaluation, not representation: a subclass implements
`_evaluate`, and everything derived from it (surface normal, area, volume,
grids, plotting) is supplied here identically for every representation.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from functools import cached_property

import numpy as np


class Surface(ABC):
    """Abstract toroidal surface, periodic in theta with period 2*pi and in
    zeta with period 2*pi/nfp.

    Attributes set by subclasses: `nfp`, `stellarator_symmetric`, `ntheta`,
    `nzeta`, `standard_toroidal_angle`.
    """

    nfp: int
    stellarator_symmetric: bool
    ntheta: int
    nzeta: int
    #: True if the surface's `zeta` parameter is the standard toroidal angle
    #: (`atan2(y, x)`), so a constant-`zeta` slice of `r` is a constant
    #: physical-toroidal-angle cross section. False for a surface built by
    #: moving another surface's points along its normal without re-solving
    #: for the standard toroidal angle (see
    #: `CoilSurface.from_uniform_offset(..., standard_toroidal_angle=False)`)
    #: -- for such a surface, `r[:, :, k]` is not a plane of constant
    #: physical toroidal angle, and cross-section plots must not assume it
    #: is.
    standard_toroidal_angle: bool

    @abstractmethod
    def _evaluate(self, theta: np.ndarray, zetal: np.ndarray) -> dict:
        """Evaluate the surface and its first derivatives.

        Parameters
        ----------
        theta : (ntheta,) array
        zetal : (nzetal,) array, the *unwrapped* toroidal angle (one full
            torus, i.e. ``nzetal = nzeta * nfp`` points spanning ``[0, 2*pi)``).

        Returns
        -------
        dict with keys 'r', 'drdtheta', 'drdzeta', each of shape
        (3, len(theta), len(zetal)) holding Cartesian (x, y, z) components.
        """
        raise NotImplementedError

    @property
    def nzetal(self) -> int:
        return self.nzeta * self.nfp

    @cached_property
    def theta(self) -> np.ndarray:
        return 2 * np.pi * np.arange(self.ntheta) / self.ntheta

    @cached_property
    def zeta(self) -> np.ndarray:
        """One field period."""
        return (2 * np.pi / self.nfp) * np.arange(self.nzeta) / self.nzeta

    @cached_property
    def zetal(self) -> np.ndarray:
        """The full torus (all field periods)."""
        return 2 * np.pi * np.arange(self.nzetal) / self.nzetal

    @cached_property
    def dtheta(self) -> float:
        return self.theta[1] - self.theta[0]

    @cached_property
    def dzeta(self) -> float:
        return self.zeta[1] - self.zeta[0]

    @cached_property
    def _evaluated(self) -> dict:
        return self._evaluate(self.theta, self.zetal)

    @property
    def r(self) -> np.ndarray:
        """(3, ntheta, nzetal) Cartesian position, all field periods."""
        return self._evaluated["r"]

    @property
    def drdtheta(self) -> np.ndarray:
        return self._evaluated["drdtheta"]

    @property
    def drdzeta(self) -> np.ndarray:
        return self._evaluated["drdzeta"]

    @cached_property
    def normal(self) -> np.ndarray:
        """(3, ntheta, nzetal) un-normalized surface normal, N = dr/dzeta x dr/dtheta."""
        drdtheta = self.drdtheta
        drdzeta = self.drdzeta
        normal = np.empty_like(drdtheta)
        normal[0] = drdzeta[1] * drdtheta[2] - drdtheta[1] * drdzeta[2]
        normal[1] = drdzeta[2] * drdtheta[0] - drdtheta[2] * drdzeta[0]
        normal[2] = drdzeta[0] * drdtheta[1] - drdtheta[0] * drdzeta[1]
        return normal

    @cached_property
    def norm_normal(self) -> np.ndarray:
        """(ntheta, nzeta) |N|, one field period."""
        n = self.normal[:, :, : self.nzeta]
        return np.sqrt(np.sum(n * n, axis=0))

    @cached_property
    def area(self) -> float:
        return self.nfp * self.dtheta * self.dzeta * np.sum(self.norm_normal)

    @cached_property
    def volume(self) -> float:
        """Enclosed volume via int (1/2) R^2 dZ dzeta (R^2 interpolated
        full -> half theta grid); `r` already spans all field periods, so no
        extra factor of nfp is needed.
        """
        r = self.r
        major_R_squared = r[0] * r[0] + r[1] * r[1]
        Z = r[2]
        interior = np.sum(
            0.5 * (major_R_squared[:-1, :] + major_R_squared[1:, :]) * (Z[1:, :] - Z[:-1, :])
        )
        wrap = np.sum(0.5 * (major_R_squared[0, :] + major_R_squared[-1, :]) * (Z[0, :] - Z[-1, :]))
        return abs((interior + wrap) * self.dzeta / 2)

    def cross_section(self, phi=None):
        """Surface cross section(s) at fixed *physical* toroidal angle(s).

        Unlike a constant-`zeta` slice, this is correct regardless of
        `standard_toroidal_angle` (ADR-025): `phi = atan2(y, x)` is computed
        from the actual Cartesian `r` grid, and each theta-line is
        interpolated (periodically in the full-torus `zetal` direction) to
        the requested physical angle(s) (ADR-029 decision 6).

        Parameters
        ----------
        phi : array-like, optional
            Physical toroidal angles in radians. Defaults to
            ``[0, 0.5, 1, 1.5] * pi / nfp`` (half a field period, four
            slices -- the legacy `regcoilPlot` set).

        Returns
        -------
        R, Z : (len(phi), ntheta) arrays.
        """
        if phi is None:
            phi = np.array([0, 0.5, 1, 1.5]) * np.pi / self.nfp
        phi = np.atleast_1d(np.asarray(phi, dtype=float))
        phi_wrapped = np.mod(phi, 2 * np.pi)

        r = self.r
        X, Y, Z = r[0], r[1], r[2]
        R = np.hypot(X, Y)
        phys_phi = np.arctan2(Y, X)  # (ntheta, nzetal), in [-pi, pi]
        phys_phi = np.mod(phys_phi, 2 * np.pi)

        ntheta = self.ntheta
        R_slices = np.empty((len(phi), ntheta))
        Z_slices = np.empty((len(phi), ntheta))
        for itheta in range(ntheta):
            ang = np.unwrap(phys_phi[itheta])
            # Extend by +/- one full turn so interpolation works across the
            # phi=0/2*pi seam regardless of where `ang` starts.
            ang_ext = np.concatenate([ang - 2 * np.pi, ang, ang + 2 * np.pi])
            R_ext = np.tile(R[itheta], 3)
            Z_ext = np.tile(Z[itheta], 3)
            order = np.argsort(ang_ext)
            R_slices[:, itheta] = np.interp(phi_wrapped, ang_ext[order], R_ext[order])
            Z_slices[:, itheta] = np.interp(phi_wrapped, ang_ext[order], Z_ext[order])
        return R_slices, Z_slices

    def plot(self, ax=None, **kwargs):
        """Minimal 3D wireframe/surface plot (matplotlib). Superseded by the
        interactive `regcoil.plot.plot_3d` (Plotly) for real use; kept as a
        lightweight fallback."""
        import matplotlib.pyplot as plt

        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(projection="3d")
        r = self.r
        kwargs.setdefault("rstride", 1)
        kwargs.setdefault("cstride", 1)
        ax.plot_surface(r[0], r[1], r[2], **kwargs)
        ax.set_box_aspect((1, 1, 1))
        return ax

    def plot_cross_section(self, other=None, phi=None, ax=None):
        """Convenience delegate. With `other=None`, single-surface overlay
        via `regcoil.plot.cross_sections_overlay(self, phi=phi, ax=ax)`
        (color by phi, returns the `ax`). With `other` given (the plasma or
        coil counterpart of `self`), the multi-subplot grid via
        `regcoil.plot.cross_sections(self, other, phi=phi)` (plasma red /
        coil blue, one subplot per phi, returns a `Figure`; `ax` is unused in
        this case since a whole grid of axes is created)."""
        from . import plot

        if other is None:
            return plot.cross_sections_overlay(self, phi=phi, ax=ax)

        from .plasma_surface import PlasmaSurface

        plasma, coil = (self, other) if isinstance(self, PlasmaSurface) else (other, self)
        return plot.cross_sections(plasma, coil, phi=phi)
