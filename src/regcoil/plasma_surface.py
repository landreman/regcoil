"""`PlasmaSurface`: the plasma boundary, plus the physics REGCOIL needs from it."""

from __future__ import annotations

import numpy as np

from . import _io
from ._constants import MU0
from ._io import bnormal_from_focus_modes, read_ascii_table, read_bnorm_file, read_focus_boundary, read_vmec_wout
from .fourier_surface import FourierSurface
from .surface import _single_target


class PlasmaSurface(FourierSurface):
    """Adds the physics attached to the plasma boundary: B_normal from the
    plasma current, and the net poloidal current (for the external-field
    boundary condition).
    """

    #: `B_normal` on the plasma surface from the plasma's own current (as
    #: opposed to the coils). Zero unless set by `from_focus` or
    #: `set_bnormal_from_bnorm_file`/`set_bnormal_from_virtual_casing`.
    Bnormal_from_plasma_current: np.ndarray
    #: `= 2*pi*B_zeta/mu0 ~ B*L/mu0`, the net poloidal current linking the
    #: plasma (for the external-field boundary condition). Set by `from_wout`;
    #: 1.0 otherwise.
    net_poloidal_current: float
    #: `= (2*pi/nfp)*bsubvmnc(0,0) ~ B*L`, multiplies the data in a bnorm
    #: file (see `set_bnormal_from_bnorm_file`). Set by `from_wout`; 1.0
    #: otherwise.
    curpol: float
    #: `|B|` on the plasma boundary. Set by `from_wout`; identically zero
    #: otherwise.
    modB: np.ndarray
    #: The volume-averaged `|B|`, or None unless set by `from_wout`.
    volume_averaged_B: float | None

    def __init__(self, *args, net_poloidal_current=1.0, curpol=1.0, **kwargs):
        super().__init__(*args, **kwargs)
        self.net_poloidal_current = net_poloidal_current
        self.curpol = curpol
        self.Bnormal_from_plasma_current = np.zeros((self.ntheta, self.nzeta))
        self.modB = np.zeros((self.ntheta, self.nzeta))
        self.volume_averaged_B = None

    @classmethod
    def from_wout(
        cls,
        wout_filename,
        ntheta=64,
        nzeta=65,
    ):
        """Build a plasma surface from a VMEC `wout` file or a simsopt `Vmec` object.

        `wout_filename` is either a path to a NetCDF `wout` file, or a
        `simsopt.mhd.Vmec`-like object (anything with a `wout` attribute).
        simsopt is not required unless you pass such an object.

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

        surf.volume_averaged_B = float(data["volavgB"])

        return surf

    @classmethod
    def from_ascii_table(cls, filename, nfp, ntheta=64, nzeta=65):
        """Read a plain ASCII-table plasma shape (legacy geometry_option_plasma=6)."""
        data = read_ascii_table(filename)
        return cls(
            data["xm"], data["xn"], data["rmnc"], data["zmns"], data["rmns"], data["zmnc"],
            nfp=nfp, ntheta=ntheta, nzeta=nzeta, stellarator_symmetric=False,
        )

    @classmethod
    def from_focus(cls, filename, ntheta=64, nzeta=65):
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

    def scale(
        self,
        *,
        length_factor=None,
        minor_radius=None,
        major_radius=None,
        volume=None,
        area=None,
        B_factor=None,
        volume_averaged_B=None,
        max_modB=None,
    ):
        """Return a new `PlasmaSurface`, scaled in size and/or in magnetic
        field strength.

        Give at most one *length* target and at most one *field* target; the
        two are independent, and either may be omitted (no scaling of that
        kind).

        Length targets (see `Surface.scale`):

        - `length_factor`: the factor itself.
        - `minor_radius`, `major_radius`, `volume`, `area`: the desired value.

        Field targets:

        - `B_factor`: the factor itself.
        - `volume_averaged_B`: the desired volume-averaged `|B|`. Only
          available if this surface has one (it is read from a VMEC `wout`
          file by `from_wout`).
        - `max_modB`: the desired maximum of `|B|` on the boundary. Only
          available if `modB` is set (again, by `from_wout`).

        `self` is never modified; a new `PlasmaSurface` is returned.

        Examples
        --------
        Scale to a desired minor radius and volume-averaged field::

            big = plasma.scale(minor_radius=1.7, volume_averaged_B=5.7)

        or give the two factors directly::

            doubled = plasma.scale(length_factor=2.0, B_factor=1.5)

        Notes
        -----
        Which quantity picks up which factor follows from its dimensions,
        with `cL` the length factor and `cB` the field factor:

        - the shape (`rmnc`, `rmns`, `zmnc`, `zmns`): `cL`
        - `net_poloidal_current` (`= 2*pi*B_zeta/mu0 ~ B*L/mu0`): `cL * cB`
        - `curpol` (`= (2*pi/nfp)*bsubvmnc(0,0) ~ B*L`): `cL * cB`
        - `modB`, `volume_averaged_B`, `Bnormal_from_plasma_current`: `cB`

        In particular the fields are *not* length-scaled: Biot-Savart gives
        `B ~ mu0*I/L = cL*cB/cL = cB`, which is why a pure size change leaves
        `|B|` alone while the current that produces it grows.

        See `Surface.scale` for the two caveats that also apply here: a
        `CoilSurface` already built from this plasma surface is not scaled
        along with it (scale it by the same `cL`, and scale a
        `from_uniform_offset` `separation` too), and a `Regcoil` copies
        `net_poloidal_current` at construction, so scale before building one.
        """
        length = self._length_scale_factor(
            length_factor=length_factor,
            minor_radius=minor_radius,
            major_radius=major_radius,
            volume=volume,
            area=area,
        )
        field = self._field_scale_factor(
            B_factor=B_factor,
            volume_averaged_B=volume_averaged_B,
            max_modB=max_modB,
        )
        return self._scaled(length, field)

    def scale_B(self, factor):
        """Return a new `PlasmaSurface` with the magnetic field strength
        multiplied by `factor`. Shorthand for `scale(B_factor=factor)`."""
        return self.scale(B_factor=factor)

    def _field_scale_factor(self, *, B_factor, volume_averaged_B, max_modB):
        """Resolve the mutually exclusive field keywords of `scale` to a
        single factor. Both field targets are linear in `B`, so the factor is
        just `target / current`."""
        key, value = _single_target(
            {
                "B_factor": B_factor,
                "volume_averaged_B": volume_averaged_B,
                "max_modB": max_modB,
            },
            "field",
        )
        if key is None:
            return 1.0
        if key == "B_factor":
            return value
        if key == "volume_averaged_B":
            if self.volume_averaged_B is None:
                raise ValueError(
                    "Cannot scale to a target volume_averaged_B: this "
                    "PlasmaSurface has no volume_averaged_B (it is set only by "
                    "from_wout). Pass B_factor instead."
                )
            current = self.volume_averaged_B
        else:
            current = float(np.max(np.abs(self.modB)))
            if current == 0:
                raise ValueError(
                    "Cannot scale to a target max_modB: this PlasmaSurface's "
                    "modB is identically zero (it is set only by from_wout). "
                    "Pass B_factor instead."
                )
        return value / current

    def _apply_length_scaling(self, factor):
        """`Surface._apply_length_scaling`, extended: on top of the shape, the
        two quantities that carry a power of length, `net_poloidal_current`
        (`~ B*L/mu0`) and `curpol` (`~ B*L`). The field quantities themselves
        (`modB`, `volume_averaged_B`, `Bnormal_from_plasma_current`) are
        length-invariant."""
        super()._apply_length_scaling(factor)
        self.net_poloidal_current *= factor
        self.curpol *= factor

    def _apply_field_scaling(self, factor):
        """`Surface._apply_field_scaling`: everything linear in `B`, which is
        the two `B*L` quantities as well as the fields themselves."""
        self.net_poloidal_current *= factor
        self.curpol *= factor
        self.modB = self.modB * factor
        self.Bnormal_from_plasma_current = self.Bnormal_from_plasma_current * factor
        if self.volume_averaged_B is not None:
            self.volume_averaged_B *= factor

    def set_bnormal_from_bnorm_file(self, filename):
        """Set `Bnormal_from_plasma_current` from a BNORM-format file.

        Requires `curpol` to already be set correctly (done automatically by
        `from_wout`).
        """
        self.Bnormal_from_plasma_current = read_bnorm_file(
            filename, self.theta, self.zeta, self.nfp, self.curpol
        )

    def set_bnormal_from_virtual_casing(self, source):
        """Set `Bnormal_from_plasma_current` from a simsopt virtual-casing result.

        An alternative to `set_bnormal_from_bnorm_file`. `source` is either the
        path to a simsopt virtual-casing NetCDF file (`vcasing*.nc`, written by
        `simsopt.mhd.VirtualCasing.save` or by the `filename` argument of
        `VirtualCasing.from_vmec`), or a `simsopt.mhd.VirtualCasing` object
        itself. Files are read directly with NetCDF, so simsopt need not be
        installed to use this.

        Unlike a BNORM file, which holds Fourier modes that can be evaluated
        anywhere, simsopt supplies `B_external_normal` on a fixed grid; it is
        interpolated (spectrally, see `_io.bnormal_from_virtual_casing`) onto
        this surface's `(theta, zeta)` grid. Both the usual `use_stellsym=True`
        layout, where the data cover half a field period, and the
        `use_stellsym=False` layout, where they cover a whole one, are accepted.

        No `curpol` scaling is involved: simsopt stores `B_external_normal` in
        Tesla, which is what REGCOIL wants, whereas BNORM writes its amplitudes
        divided by `curpol`. The two agree in sign and magnitude, so the choice
        of source is invisible downstream (see simsopt's
        `tests/mhd/test_virtual_casing.py::test_bnorm_benchmark`).

        The virtual-casing calculation must have been run on *this* plasma
        boundary; only `nfp` can be checked here.
        """
        data = _io.virtual_casing_data(source)
        self.Bnormal_from_plasma_current = _io.bnormal_from_virtual_casing(
            data, self.theta, self.zeta, self.nfp
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
