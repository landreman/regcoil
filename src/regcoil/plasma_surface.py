"""`PlasmaSurface`: the plasma boundary, plus the physics REGCOIL needs from it."""

from __future__ import annotations

import numpy as np

from ._constants import MU0
from ._io import bnormal_from_focus_modes, read_ascii_table, read_bnorm_file, read_focus_boundary, read_vmec_wout
from .fourier_surface import FourierSurface


class PlasmaSurface(FourierSurface):
    """Adds the physics attached to the plasma boundary: B_normal from the
    plasma current, and the net poloidal current (for the external-field
    boundary condition).
    """

    def __init__(self, *args, net_poloidal_current=1.0, curpol=1.0, **kwargs):
        super().__init__(*args, **kwargs)
        self.net_poloidal_current = net_poloidal_current
        self.curpol = curpol
        self.Bnormal_from_plasma_current = np.zeros((self.ntheta, self.nzeta))

    @classmethod
    def from_wout(
        cls,
        wout_filename,
        ntheta=64,
        nzeta=64,
        mesh="full",
        straight_field_line=False,
    ):
        """Build a plasma surface from a VMEC `wout` file.

        `mesh="full"` uses the outermost point of VMEC's full radial mesh
        (legacy geometry_option_plasma=2); `mesh="half"` averages the two
        outermost full-mesh points to land on the half mesh
        (geometry_option_plasma=3). `straight_field_line=True` instead
        transforms to the straight-field-line poloidal angle
        (geometry_option_plasma=4; not implemented for non-stellarator-
        symmetric equilibria, matching the legacy limitation).
        """
        data = read_vmec_wout(wout_filename)
        nfp = data["nfp"]
        lasym = data["lasym"]

        if straight_field_line:
            # Not implemented: the legacy root-solve (regcoil_fzero bracketing
            # VMEC's theta against a fixed +/-0.3 window) is fragile even in
            # the reference Fortran -- it fails to bracket a root ("no sign
            # change in residual") for some equilibria/grids -- and there is
            # no validated reference to port against.
            raise NotImplementedError(
                "from_wout(straight_field_line=True) is not implemented (legacy "
                "geometry_option_plasma=4). The legacy root-solve is not robust "
                "enough to port with confidence; use mesh='full' or mesh='half' "
                "instead, or supply your own straight-field-line coefficients "
                "directly via PlasmaSurface(...)."
            )
        else:
            if mesh == "full":
                weight1, weight2 = 0.0, 1.0
            elif mesh == "half":
                weight1, weight2 = 0.5, 0.5
            else:
                raise ValueError(f"mesh must be 'full' or 'half', got {mesh!r}")

            rmnc = data["rmnc"][-2] * weight1 + data["rmnc"][-1] * weight2
            zmns = data["zmns"][-2] * weight1 + data["zmns"][-1] * weight2
            rmns = zmnc = None
            if lasym:
                rmns = data["rmns"][-2] * weight1 + data["rmns"][-1] * weight2
                zmnc = data["zmnc"][-2] * weight1 + data["zmnc"][-1] * weight2

            surf = cls(
                data["xm"], data["xn"], rmnc, zmns, rmns, zmnc,
                nfp=nfp, ntheta=ntheta, nzeta=nzeta, stellarator_symmetric=not lasym,
            )

        # VMEC stores the toroidal Boozer covariant component B_zeta ("bvco")
        # on the half mesh; extrapolate to the boundary the same way the
        # legacy code does (1.5*last - 0.5*second-to-last).
        surf.net_poloidal_current = (
            2 * np.pi / MU0 * (1.5 * data["bvco"][-1] - 0.5 * data["bvco"][-2])
        )
        # curpol multiplies the data in a bnorm file (see set_bnormal_from_bnorm_file).
        surf.curpol = (2 * np.pi / nfp) * (
            1.5 * data["bsubvmnc"][-1, 0] - 0.5 * data["bsubvmnc"][-2, 0]
        )
        return surf

    @classmethod
    def from_ascii_table(cls, filename, nfp, ntheta=64, nzeta=64):
        """Read a plain ASCII-table plasma shape (legacy geometry_option_plasma=6)."""
        data = read_ascii_table(filename)
        return cls(
            data["xm"], data["xn"], data["rmnc"], data["zmns"], data["rmns"], data["zmnc"],
            nfp=nfp, ntheta=ntheta, nzeta=nzeta, stellarator_symmetric=False,
        )

    @classmethod
    def from_focus(cls, filename, ntheta=64, nzeta=64):
        """Read a FOCUS-format plasma boundary; also sets `Bnormal_from_plasma_current`
        from the file's Bn Fourier modes, if present.
        """
        data = read_focus_boundary(filename)
        surf = cls(
            data["xm"], data["xn"], data["rmnc"], data["zmns"], data["rmns"], data["zmnc"],
            nfp=data["nfp"], ntheta=ntheta, nzeta=nzeta, stellarator_symmetric=False,
        )
        if data["nbf"] > 0:
            surf.Bnormal_from_plasma_current = bnormal_from_focus_modes(
                data["bfm"], data["bfn"], data["bfc"], data["bfs"], surf.theta, surf.zeta
            )
        return surf

    def set_bnormal_from_bnorm_file(self, filename):
        """Set `Bnormal_from_plasma_current` from a BNORM-format file.

        Requires `curpol` to already be set correctly (done automatically by
        `from_wout`).
        """
        self.Bnormal_from_plasma_current = read_bnorm_file(
            filename, self.theta, self.zeta, self.nfp, self.curpol
        )
