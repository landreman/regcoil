"""`CoilSurface`: the coil winding surface."""

from __future__ import annotations

import numpy as np

from ._io import read_nescin
from .fourier_surface import FourierSurface


class CoilSurface(FourierSurface):
    """Nearly bare: coil-side Fourier filtering, plus
    `from_uniform_offset`, the one constructor that calls Fortran.
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
        ntheta_transform=None, nzeta_transform=None, tol=1e-10,
    ):
        """A surface offset uniformly outward from `plasma` by `separation`
        meters along its normal. `ntheta_transform`/`nzeta_transform` (default
        `ntheta`/`nzeta`) set the grid the offset surface is sampled on before
        being Fourier-transformed to `mpol`/`ntor` modes; each grid point is a
        Fortran root solve (`regcoil_uniform_offset_surface`, Phase 7).
        """
        from . import _core

        ntheta_transform = ntheta if ntheta_transform is None else ntheta_transform
        nzeta_transform = nzeta if nzeta_transform is None else nzeta_transform

        xm_out, xn_out, rmnc_out, rmns_out, zmnc_out, zmns_out = _core.uniform_offset_surface(
            plasma.xm, plasma.xn, plasma.rmnc, plasma.rmns, plasma.zmnc, plasma.zmns,
            not plasma.stellarator_symmetric, plasma.nfp,
            float(separation), int(mpol), int(ntor),
            int(ntheta_transform), int(nzeta_transform), float(tol),
        )
        return cls(
            xm_out, xn_out, rmnc_out, zmns_out, rmns_out, zmnc_out,
            nfp=plasma.nfp, ntheta=ntheta, nzeta=nzeta,
        )

    def filter_modes(self, mpol_filter=None, ntor_filter=None):
        """Zero out modes with |m| > mpol_filter or |n| > ntor_filter*nfp, in place."""
        mpol_limit = np.inf if mpol_filter is None else mpol_filter
        ntor_limit = np.inf if ntor_filter is None else ntor_filter * self.nfp
        mask = (np.abs(self.xm) > mpol_limit) | (np.abs(self.xn) > ntor_limit)
        self.rmnc[mask] = 0
        self.rmns[mask] = 0
        self.zmnc[mask] = 0
        self.zmns[mask] = 0
