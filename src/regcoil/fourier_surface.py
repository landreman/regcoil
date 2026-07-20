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

from .surface import Surface


class FourierSurface(Surface):
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

    def _evaluate(self, theta, zetal):
        # The double-angle expansion turns the (mnmax, ntheta, nzetal) sum
        # into a handful of (ntheta, mnmax) @ (mnmax, nzetal) gemms.
        m = self.xm.astype(np.float64)[:, None]
        n = self.xn.astype(np.float64)[:, None]
        cos_mtheta = np.cos(np.outer(self.xm, theta))  # (mnmax, ntheta)
        sin_mtheta = np.sin(np.outer(self.xm, theta))
        cos_nzeta = np.cos(np.outer(self.xn, zetal))  # (mnmax, nzetal)
        sin_nzeta = np.sin(np.outer(self.xn, zetal))

        rmnc = self.rmnc[:, None]
        rmns = self.rmns[:, None]
        zmnc = self.zmnc[:, None]
        zmns = self.zmns[:, None]

        # R = P1^T @ cos_nzeta + P2^T @ sin_nzeta; Z likewise with Q1, Q2.
        P1 = rmnc * cos_mtheta + rmns * sin_mtheta
        P2 = rmnc * sin_mtheta - rmns * cos_mtheta
        Q1 = zmns * sin_mtheta + zmnc * cos_mtheta
        Q2 = zmnc * sin_mtheta - zmns * cos_mtheta

        dP1dtheta = -m * P2
        dP2dtheta = m * P1
        dQ1dtheta = -m * Q2
        dQ2dtheta = m * Q1

        dcos_nzeta = -n * sin_nzeta
        dsin_nzeta = n * cos_nzeta

        R = P1.T @ cos_nzeta + P2.T @ sin_nzeta
        Z = Q1.T @ cos_nzeta + Q2.T @ sin_nzeta
        dRdtheta = dP1dtheta.T @ cos_nzeta + dP2dtheta.T @ sin_nzeta
        dZdtheta = dQ1dtheta.T @ cos_nzeta + dQ2dtheta.T @ sin_nzeta
        dRdzeta = P1.T @ dcos_nzeta + P2.T @ dsin_nzeta
        dZdzeta = Q1.T @ dcos_nzeta + Q2.T @ dsin_nzeta

        # Toroidal angle shift nu (same mode structure as Z): phi = zeta + nu.
        # dphidtheta = dnudtheta, dphidzeta = 1 + dnudzeta. The common case
        # (nu identically zero -> phi = zeta) keeps the original fast path.
        if self.standard_toroidal_angle:
            nu = 0.0
            dnudtheta = 0.0
            dnudzeta = 0.0
            phi = zetal[None, :]
        else:
            numns = self.numns[:, None]
            numnc = self.numnc[:, None]
            N1 = numns * sin_mtheta + numnc * cos_mtheta
            N2 = numnc * sin_mtheta - numns * cos_mtheta
            dN1dtheta = -m * N2
            dN2dtheta = m * N1
            nu = N1.T @ cos_nzeta + N2.T @ sin_nzeta
            dnudtheta = dN1dtheta.T @ cos_nzeta + dN2dtheta.T @ sin_nzeta
            dnudzeta = N1.T @ dcos_nzeta + N2.T @ dsin_nzeta
            phi = zetal[None, :] + nu

        cosphi = np.cos(phi)
        sinphi = np.sin(phi)

        X = R * cosphi
        Y = R * sinphi
        dXdtheta = dRdtheta * cosphi - R * sinphi * dnudtheta
        dYdtheta = dRdtheta * sinphi + R * cosphi * dnudtheta
        dXdzeta = dRdzeta * cosphi - R * sinphi * (1.0 + dnudzeta)
        dYdzeta = dRdzeta * sinphi + R * cosphi * (1.0 + dnudzeta)

        r = np.stack([X, Y, Z], axis=0)
        drdtheta = np.stack([dXdtheta, dYdtheta, dZdtheta], axis=0)
        drdzeta = np.stack([dXdzeta, dYdzeta, dZdzeta], axis=0)
        return {"r": r, "drdtheta": drdtheta, "drdzeta": drdzeta}

    def evaluate_at(self, theta_pts, zeta_pts):
        """Evaluate the surface (and its first theta/zeta derivatives) at
        arbitrary *paired* `(theta, zeta)` points, rather than the
        tensor-product grid `_evaluate`/`r` use. Used where points come from
        a contour (`regcoil.cut`) rather than a regular grid.

        Returns a dict with keys 'r', 'drdtheta', 'drdzeta', each
        `(3, len(theta_pts))`, matching `_evaluate`'s per-point contract.
        """
        theta_pts = np.asarray(theta_pts, dtype=np.float64)
        zeta_pts = np.asarray(zeta_pts, dtype=np.float64)
        m = self.xm[:, None].astype(np.float64)
        n = self.xn[:, None].astype(np.float64)
        angle = self.xm[:, None] * theta_pts[None, :] - self.xn[:, None] * zeta_pts[None, :]
        cos_a = np.cos(angle)
        sin_a = np.sin(angle)

        rmnc, rmns = self.rmnc[:, None], self.rmns[:, None]
        zmnc, zmns = self.zmnc[:, None], self.zmns[:, None]

        R = np.sum(rmnc * cos_a + rmns * sin_a, axis=0)
        Z = np.sum(zmns * sin_a + zmnc * cos_a, axis=0)
        R_term = -rmnc * sin_a + rmns * cos_a
        Z_term = zmns * cos_a - zmnc * sin_a
        dRdtheta = np.sum(m * R_term, axis=0)
        dRdzeta = np.sum(-n * R_term, axis=0)
        dZdtheta = np.sum(m * Z_term, axis=0)
        dZdzeta = np.sum(-n * Z_term, axis=0)

        if self.standard_toroidal_angle:
            phi = zeta_pts
            dphidtheta = np.zeros_like(theta_pts)
            dphidzeta = np.ones_like(zeta_pts)
        else:
            numns, numnc = self.numns[:, None], self.numnc[:, None]
            nu = np.sum(numns * sin_a + numnc * cos_a, axis=0)
            nu_term = numns * cos_a - numnc * sin_a
            dnudtheta = np.sum(m * nu_term, axis=0)
            dnudzeta = np.sum(-n * nu_term, axis=0)
            phi = zeta_pts + nu
            dphidtheta = dnudtheta
            dphidzeta = 1.0 + dnudzeta

        cosphi, sinphi = np.cos(phi), np.sin(phi)
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
