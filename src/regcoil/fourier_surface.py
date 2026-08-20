"""`FourierSurface`: a toroidal surface given by a double Fourier series.

    R(theta, zeta)  = sum_mn rmnc_mn*cos(m*theta - n*zeta) + rmns_mn*sin(m*theta - n*zeta)
    Z(theta, zeta)  = sum_mn zmns_mn*sin(m*theta - n*zeta) + zmnc_mn*cos(m*theta - n*zeta)
    nu(theta, zeta) = sum_mn numns_mn*sin(m*theta - n*zeta) + numnc_mn*cos(m*theta - n*zeta)

The Cartesian point is placed at physical toroidal angle ``phi = zeta + nu``:

    x = R*cos(zeta + nu),  y = R*sin(zeta + nu),  z = Z.

``nu`` is the deviation of the standard (physical) toroidal angle from the
surface parameter ``zeta``. When all ``nu`` modes are zero (the usual case,
and the default) ``phi = zeta``, ``standard_toroidal_angle`` is True, and this
reduces to an ordinary VMEC-style Fourier surface. A nonzero ``nu`` lets a
surface whose points are *not* uniformly spaced in the standard toroidal angle
-- e.g. the normal-offset coil surface from
`CoilSurface.from_uniform_offset(..., standard_toroidal_angle=False)` -- still
be represented (and smoothed / mode-filtered) as a Fourier surface. Under
stellarator symmetry ``nu`` has sine parity (``numnc = 0``), like ``Z``'s
``zmnc``.

Conventions (asserted below; see API.md): `xn` already includes the factor of
`nfp` (VMEC convention), and the angle is `m*theta - n*zeta`.
"""

from __future__ import annotations

import numpy as np

from ._transforms import evaluate_modes, evaluate_modes_by_column
from .surface import Surface


class FourierSurface(Surface):
    """A toroidal surface given by a double Fourier series in `theta`/`zeta`
    (see the module docstring for the exact `R`/`Z`/`nu` expansion and the
    `standard_toroidal_angle` convention)."""

    def __init__(
        self,
        xm,
        xn,
        rmnc,
        zmns,
        rmns=None,
        zmnc=None,
        *,
        nfp,
        ntheta=64,
        nzeta=64,
        stellarator_symmetric=None,
        standard_toroidal_angle=True,
        numns=None,
        numnc=None,
    ):
        xm = np.asarray(xm, dtype=np.int64)
        xn = np.asarray(xn, dtype=np.int64)
        rmnc = np.asarray(rmnc, dtype=np.float64)
        zmns = np.asarray(zmns, dtype=np.float64)
        mnmax = xm.shape[0]
        for name, arr in (("xn", xn), ("rmnc", rmnc), ("zmns", zmns)):
            if arr.shape != (mnmax,):
                raise ValueError(f"{name} must have shape (mnmax,) = ({mnmax},), got {arr.shape}")

        rmns = np.zeros(mnmax) if rmns is None else np.asarray(rmns, dtype=np.float64)
        zmnc = np.zeros(mnmax) if zmnc is None else np.asarray(zmnc, dtype=np.float64)
        numns = np.zeros(mnmax) if numns is None else np.asarray(numns, dtype=np.float64)
        numnc = np.zeros(mnmax) if numnc is None else np.asarray(numnc, dtype=np.float64)
        for name, arr in (("rmns", rmns), ("zmnc", zmnc), ("numns", numns), ("numnc", numnc)):
            if arr.shape != (mnmax,):
                raise ValueError(f"{name} must have shape (mnmax,) = ({mnmax},), got {arr.shape}")

        nfp = int(nfp)
        if nfp < 1:
            raise ValueError(f"nfp must be a positive integer, got {nfp}")
        if mnmax > 0 and np.any(xn % nfp != 0):
            raise ValueError("xn must be an integer multiple of nfp (VMEC convention)")

        if stellarator_symmetric is None:
            stellarator_symmetric = bool(
                np.all(rmns == 0) and np.all(zmnc == 0) and np.all(numnc == 0)
            )

        self.nfp = nfp
        self.ntheta = int(ntheta)
        self.nzeta = int(nzeta)
        self.stellarator_symmetric = bool(stellarator_symmetric)
        self.standard_toroidal_angle = bool(standard_toroidal_angle)
        self.mnmax = mnmax
        self.xm = xm
        self.xn = xn
        self.rmnc = rmnc
        self.rmns = rmns
        self.zmnc = zmnc
        self.zmns = zmns
        self.numns = numns
        self.numnc = numnc

    def _amplitudes(self):
        """The Fourier amplitudes of everything an evaluation needs -- `R`,
        `Z`, and (unless `zeta` is the standard toroidal angle) `nu`, each
        followed by its `theta`- and `zeta`-derivative -- as the rows of one
        `(3*nfields, mnmax)` cosine array and one sine array.

        Every one of those quantities is a series
        `sum_mn cmn*cos(m*theta - n*zeta) + smn*sin(...)` in the same modes,
        differing only in its amplitudes: differentiating gives the `d/dtheta`
        pair `(m*smn, -m*cmn)` and the `d/dzeta` pair `(-n*smn, n*cmn)`. So
        stacking them lets one pass over the grid produce them all.
        """
        m = self.xm.astype(np.float64)
        n = self.xn.astype(np.float64)
        fields = [(self.rmnc, self.rmns), (self.zmnc, self.zmns)]
        if not self.standard_toroidal_angle:
            fields.append((self.numnc, self.numns))
        cmn = np.array([row for fc, fs in fields for row in (fc, m * fs, -n * fs)])
        smn = np.array([row for fc, fs in fields for row in (fs, -m * fc, n * fc)])
        return cmn, smn

    def _cartesian(self, values, zeta):
        """Assemble the `r`/`drdtheta`/`drdzeta` dict from the stacked
        `R`/`Z`/`nu` values and derivatives returned for `_amplitudes`, at
        points whose toroidal parameter is `zeta` (broadcast against them).
        """
        R, dRdtheta, dRdzeta = values[0], values[1], values[2]
        Z, dZdtheta, dZdzeta = values[3], values[4], values[5]

        # Toroidal angle shift nu (same mode structure as Z): phi = zeta + nu,
        # so dphidtheta = dnudtheta and dphidzeta = 1 + dnudzeta. The common
        # case (nu identically zero -> phi = zeta) skips the nu series
        # entirely, in `_amplitudes`.
        if self.standard_toroidal_angle:
            phi = zeta
            dphidtheta = 0.0
            dphidzeta = 1.0
        else:
            nu, dnudtheta, dnudzeta = values[6], values[7], values[8]
            phi = zeta + nu
            dphidtheta = dnudtheta
            dphidzeta = 1.0 + dnudzeta

        cosphi = np.cos(phi)
        sinphi = np.sin(phi)

        X = R * cosphi
        Y = R * sinphi
        dXdtheta = dRdtheta * cosphi - R * sinphi * dphidtheta
        dYdtheta = dRdtheta * sinphi + R * cosphi * dphidtheta
        dXdzeta = dRdzeta * cosphi - R * sinphi * dphidzeta
        dYdzeta = dRdzeta * sinphi + R * cosphi * dphidzeta

        r = np.stack([X, Y, Z], axis=0)
        drdtheta = np.stack([dXdtheta, dYdtheta, dZdtheta], axis=0)
        drdzeta = np.stack([dXdzeta, dYdzeta, dZdzeta], axis=0)
        return {"r": r, "drdtheta": drdtheta, "drdzeta": drdzeta}

    def _evaluate(self, theta, zetal):
        # The angle `m*theta - n*zeta` separates, so the (mnmax, ntheta,
        # nzetal) sum never has to be built; see `_transforms`.
        cmn, smn = self._amplitudes()
        values = evaluate_modes(cmn, smn, self.xm, self.xn, theta, zetal)
        return self._cartesian(values, zetal[None, :])

    def _evaluate_columns(self, theta, zeta):
        """`Surface._evaluate_columns`, overridden: a shared `zeta` grid is
        still enough for the angle to separate, so this costs little more than
        `_evaluate` rather than as much as `evaluate_at`.
        """
        cmn, smn = self._amplitudes()
        values = evaluate_modes_by_column(cmn, smn, self.xm, self.xn, theta, zeta)
        return self._cartesian(values, zeta[None, :])

    def evaluate_at(self, theta_pts, zeta_pts):
        """Evaluate the surface (and its first theta/zeta derivatives) at
        arbitrary *paired* `(theta, zeta)` points, rather than the
        tensor-product grid `_evaluate`/`r` use. Used where points come from
        a contour (`regcoil.cut`) rather than a regular grid.

        Returns a dict with keys 'r', 'drdtheta', 'drdzeta', each
        `(3, len(theta_pts))`, matching `_evaluate`'s per-point contract.

        Points this general leave nothing to separate, so this is the one
        evaluator that does pay for a `(mnmax, npoints)` array of angles;
        prefer `_evaluate` or `_evaluate_columns` when the points are
        structured.
        """
        theta_pts = np.asarray(theta_pts, dtype=np.float64)
        zeta_pts = np.asarray(zeta_pts, dtype=np.float64)
        angle = self.xm[:, None] * theta_pts[None, :] - self.xn[:, None] * zeta_pts[None, :]
        cos_a = np.cos(angle)
        sin_a = np.sin(angle, out=angle)  # `angle` is dead after this

        cmn, smn = self._amplitudes()
        values = cmn @ cos_a + smn @ sin_a
        return self._cartesian(values, zeta_pts)

    def _apply_length_scaling(self, factor):
        """`Surface._apply_length_scaling` for a Fourier surface: `R` and `Z`
        carry one power of length each, while `nu` is an angle and so is
        scale-invariant. Multiplication is out-of-place because two amplitude
        arrays may legitimately be the same object (e.g. the same array passed
        for `rmns` and `zmnc`), which an in-place `*=` would scale twice.
        """
        self.rmnc = self.rmnc * factor
        self.rmns = self.rmns * factor
        self.zmnc = self.zmnc * factor
        self.zmns = self.zmns * factor

    @classmethod
    def circular_torus(cls, R0, a, nfp, ntheta=64, nzeta=64):
        """A plain circular-cross-section torus of major radius R0, minor radius a."""
        xm = np.array([0, 1])
        xn = np.array([0, 0])
        rmnc = np.array([R0, a], dtype=np.float64)
        zmns = np.array([0.0, a], dtype=np.float64)
        return cls(
            xm, xn, rmnc, zmns,
            nfp=nfp, ntheta=ntheta, nzeta=nzeta, stellarator_symmetric=True,
        )
