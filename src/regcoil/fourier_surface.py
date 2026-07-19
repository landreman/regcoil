"""`FourierSurface`: a toroidal surface given by a double Fourier series.

    R(theta, zeta) = sum_mn rmnc_mn*cos(m*theta - n*zeta) + rmns_mn*sin(m*theta - n*zeta)
    Z(theta, zeta) = sum_mn zmns_mn*sin(m*theta - n*zeta) + zmnc_mn*cos(m*theta - n*zeta)

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
        for name, arr in (("rmns", rmns), ("zmnc", zmnc)):
            if arr.shape != (mnmax,):
                raise ValueError(f"{name} must have shape (mnmax,) = ({mnmax},), got {arr.shape}")

        nfp = int(nfp)
        if nfp < 1:
            raise ValueError(f"nfp must be a positive integer, got {nfp}")
        if mnmax > 0 and np.any(xn % nfp != 0):
            raise ValueError("xn must be an integer multiple of nfp (VMEC convention)")

        if stellarator_symmetric is None:
            stellarator_symmetric = bool(np.all(rmns == 0) and np.all(zmnc == 0))

        self.nfp = nfp
        self.ntheta = int(ntheta)
        self.nzeta = int(nzeta)
        self.stellarator_symmetric = bool(stellarator_symmetric)
        self.mnmax = mnmax
        self.xm = xm
        self.xn = xn
        self.rmnc = rmnc
        self.rmns = rmns
        self.zmnc = zmnc
        self.zmns = zmns

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

        cosphi = np.cos(zetal)[None, :]
        sinphi = np.sin(zetal)[None, :]

        X = R * cosphi
        Y = R * sinphi
        dXdtheta = dRdtheta * cosphi
        dYdtheta = dRdtheta * sinphi
        dXdzeta = dRdzeta * cosphi - R * sinphi
        dYdzeta = dRdzeta * sinphi + R * cosphi

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
