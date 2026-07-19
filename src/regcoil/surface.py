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
    `nzeta`.
    """

    nfp: int
    stellarator_symmetric: bool
    ntheta: int
    nzeta: int

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

    def plot(self, ax=None, **kwargs):
        """Minimal 3D wireframe/surface plot (matplotlib)."""
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
