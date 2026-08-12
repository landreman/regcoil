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
    #: The `ThetaMap` this surface was reparameterized by, or None if its
    #: poloidal angle is whatever its constructor produced. Set by
    #: `reparameterize_theta` and by
    #: `CoilSurface.from_uniform_offset(theta_reparameterization=...)`.
    theta_map = None

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

    def evaluate_at(self, theta_pts, zeta_pts) -> dict:
        """Evaluate the surface at arbitrary *paired* `(theta, zeta)` points,
        rather than the tensor-product grid `_evaluate` uses.

        Returns a dict with keys 'r', 'drdtheta', 'drdzeta', each
        `(3, len(theta_pts))`. Needed wherever points do not form a grid: a
        contour (`regcoil.cut`), or a theta-reparameterized sampling, where
        `theta` depends on `zeta`. There is no representation-independent way
        to supply this from `_evaluate` alone, so subclasses must override it
        (`FourierSurface` does).
        """
        raise NotImplementedError(
            f"{type(self).__name__} does not implement `evaluate_at`, which is "
            f"required for evaluation at paired (theta, zeta) points."
        )

    def reparameterize_theta(self, scheme, *, mpol=None, ntor=None, ntheta=None, nzeta=None):
        """Return a `FourierSurface` describing the *same physical surface* as
        `self`, with its poloidal angle reparameterized according to `scheme`
        (see `regcoil.reparameterize`).

        `scheme` is a `UniformArclength` or `CurvatureWeighted` instance, or
        the shorthand string `"uniform_arclength"` / `"curvature"`.

        `ntheta`/`nzeta` (default `self.ntheta`/`self.nzeta`) set the grid the
        reparameterized surface is resampled on before being transformed to
        `mpol`/`ntor` modes (defaulting to that grid's Nyquist limit).

        Note that this **resamples and refits**, so the result carries one more
        Fourier truncation than `self` did -- and reparameterizing generally
        *widens* the spectrum, because composing a band-limited `R(theta)` with
        a non-affine `g` is not band-limited. Expect to need a higher `mpol`
        than the original surface used, especially for `CurvatureWeighted`.
        Reparameterizing a VMEC plasma surface is therefore usually a net loss:
        VMEC's own angle is already chosen for good spectral width. The win is
        on a surface whose angle was inherited rather than chosen -- above all
        the uniform-offset coil surface, where
        `CoilSurface.from_uniform_offset(..., theta_reparameterization=...)`
        should be preferred anyway, since it applies the map at the sampling
        stage and so fits only once, in the good angle.
        """
        from .coil_surface import _fourier_transform_offset_surface
        from .fourier_surface import FourierSurface
        from .reparameterize import theta_map

        ntheta = self.ntheta if ntheta is None else int(ntheta)
        nzeta = self.nzeta if nzeta is None else int(nzeta)
        mpol = ntheta // 2 if mpol is None else int(mpol)
        ntor = nzeta // 2 if ntor is None else int(ntor)

        def curve(theta, zeta):
            return self._evaluate(theta, zeta)["r"]

        tmap = theta_map(
            curve, scheme, nfp=self.nfp, ntheta=ntheta, nzeta=nzeta,
            stellarator_symmetric=self.stellarator_symmetric,
        )

        theta_grid = 2 * np.pi * np.arange(ntheta) / ntheta
        zeta_grid = (2 * np.pi / self.nfp) * np.arange(nzeta) / nzeta
        theta_old = tmap(theta_grid, zeta_grid)
        zeta_2d = np.broadcast_to(zeta_grid, theta_old.shape)

        r = self.evaluate_at(np.ravel(theta_old), np.ravel(zeta_2d))["r"]
        r = r.reshape((3,) + theta_old.shape)
        major_r = np.hypot(r[0], r[1])
        nu_val = np.arctan2(r[1], r[0]) - zeta_grid[None, :]
        nu_val = np.mod(nu_val + np.pi, 2 * np.pi) - np.pi

        lasym = not self.stellarator_symmetric
        (
            xm, xn, rmnc, rmns, zmnc, zmns, numns, numnc,
        ) = _fourier_transform_offset_surface(
            theta_grid, zeta_grid, major_r, r[2], nu_val, self.nfp, mpol, ntor, lasym
        )
        if self.standard_toroidal_angle:
            # A theta-map does not touch zeta or phi, so phi = zeta survives it.
            numns = numnc = None

        surface = FourierSurface(
            xm, xn, rmnc, zmns, rmns, zmnc,
            nfp=self.nfp, ntheta=ntheta, nzeta=nzeta,
            stellarator_symmetric=self.stellarator_symmetric,
            standard_toroidal_angle=self.standard_toroidal_angle,
            numns=numns, numnc=numnc,
        )
        surface.theta_map = tmap
        return surface

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
        """(ntheta, nzeta) `|N|`, one field period."""
        n = self.normal[:, :, : self.nzeta]
        return np.sqrt(np.sum(n * n, axis=0))

    @cached_property
    def area(self) -> float:
        return self.nfp * self.dtheta * self.dzeta * np.sum(self.norm_normal)

    @cached_property
    def volume(self) -> float:
        """Enclosed volume via the coordinate-free surface integral
        `V = |int Z (dr/dzeta x dr/dtheta)_z dtheta dzeta|` (the divergence
        theorem applied to the field `(0, 0, Z)`, whose divergence is 1).

        This form is valid for *any* parameterization, including a surface
        whose `zeta` is not the standard toroidal angle (nonzero `nu`); the
        simpler `int (1/2) R^2 (dZ/dtheta) dtheta dzeta` is only correct when
        `zeta == phi`, and the two agree to machine precision in that case.
        The integrand is smooth and periodic in both theta and zeta, so the
        plain periodic sum (the periodic trapezoid rule) converges spectrally,
        and the analytic derivatives from `drdtheta`/`drdzeta` are used rather
        than a finite difference of the grid. `r` already spans all field
        periods, so no extra factor of nfp is needed.
        """
        Z = self.r[2]
        drdtheta = self.drdtheta
        drdzeta = self.drdzeta
        # z-component of the (un-normalized) surface normal dr/dzeta x dr/dtheta:
        normal_z = drdzeta[0] * drdtheta[1] - drdzeta[1] * drdtheta[0]
        return abs(np.sum(Z * normal_z) * self.dtheta * self.dzeta)

    @cached_property
    def mean_cross_sectional_area(self) -> float:
        r"""Mean cross-sectional area, averaged over the physical toroidal
        angle:

        .. math::
            \overline{A} = \frac{1}{2\pi} \int_0^{2\pi} \left(
                \int_{S_\phi} dR\, dZ \right) d\phi,

        matching simsopt's `Surface.mean_cross_sectional_area`. The change of
        variables from the physical angle `phi` to the surface's own `zeta`
        parameter turns this into an average over the `(theta, zeta)` grid, so
        it is correct even when `zeta` is not the standard toroidal angle.
        """
        r = self.r
        drdtheta = self.drdtheta
        drdzeta = self.drdzeta
        x, y = r[0], r[1]
        x2y2 = x * x + y * y
        # d(phi)/d(zeta) and d(phi)/d(theta) with phi = atan2(y, x):
        dphi_dzeta = (x * drdzeta[1] - y * drdzeta[0]) / x2y2
        dphi_dtheta = (x * drdtheta[1] - y * drdtheta[0]) / x2y2
        # R * dZ/d(theta) * |d(phi,theta)/d(zeta,theta)| after change of vars:
        integrand = np.sqrt(x2y2) * (drdtheta[2] * dphi_dzeta - drdzeta[2] * dphi_dtheta)
        return 2 * np.pi * abs(np.mean(integrand))

    @cached_property
    def minor_radius(self) -> float:
        r"""Minor radius, :math:`R_\text{minor} = \sqrt{\overline{A} / \pi}`,
        where :math:`\overline{A}` is `mean_cross_sectional_area`. Same
        definition as simsopt and VMEC."""
        return np.sqrt(self.mean_cross_sectional_area / np.pi)

    @cached_property
    def major_radius(self) -> float:
        r"""Major radius,
        :math:`R_\text{major} = V / (2 \pi^2 R_\text{minor}^2)`, where `V` is
        the enclosed `volume`. Same definition as simsopt and VMEC."""
        return abs(self.volume) / (2 * np.pi**2 * self.minor_radius**2)

    @cached_property
    def aspect_ratio(self) -> float:
        r"""Aspect ratio,
        :math:`R_\text{major} / R_\text{minor}`, using the VMEC definition.
        Same definition as simsopt."""
        return self.major_radius / self.minor_radius

    def cross_section(self, phi=None):
        """Surface cross section(s) at fixed *physical* toroidal angle(s).

        Unlike a constant-`zeta` slice, this is correct regardless of
        `standard_toroidal_angle`: `phi = atan2(y, x)` is computed
        from the actual Cartesian `r` grid, and each theta-line is
        interpolated (periodically in the full-torus `zetal` direction) to
        the requested physical angle(s).

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
