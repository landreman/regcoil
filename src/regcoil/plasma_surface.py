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
        self.modB = np.zeros((self.ntheta, self.nzeta))

    @classmethod
    def from_wout(
        cls,
        wout_filename,
        ntheta=64,
        nzeta=64,
    ):
        """Build a plasma surface from a VMEC `wout` file.

        Uses the outermost point of VMEC's full radial mesh
        (legacy geometry_option_plasma=2).
        """
        data = read_vmec_wout(wout_filename)
        nfp = data["nfp"]
        lasym = data["lasym"]

        rmnc = data["rmnc"][-1]
        zmns = data["zmns"][-1]
        # Extrapolate from half grid to boundary
        bmnc = 1.5 * data["bmnc"][-1] - 0.5 * data["bmnc"][-2]
        rmns = zmnc = bmns = None
        if lasym:
            rmns = data["rmns"][-1]
            zmnc = data["zmnc"][-1]
            bmns = 1.5 * data["bmns"][-1] - 0.5 * data["bmns"][-2]

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
        # |B| on the boundary: same cosine/sine Fourier series as R, using the
        # (possibly denser) Nyquist mode set that VMEC stores bmnc/bmns on.
        if bmns is None:
            bmns = np.zeros_like(bmnc)
        angle = (
            data["xm_nyq"][:, None, None] * surf.theta[None, :, None]
            - data["xn_nyq"][:, None, None] * surf.zeta[None, None, :]
        )
        surf.modB = (
            bmnc[:, None, None] * np.cos(angle) + bmns[:, None, None] * np.sin(angle)
        ).sum(axis=0)

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

    def save(self, path):
        """Save this plasma surface to `path`."""
        from . import _serialize
        _serialize.save(path, plasma=self)

    @classmethod
    def load(cls, path):
        """Load a `PlasmaSurface` from `path`."""
        from . import _serialize
        return _serialize.load_plasma(path)
