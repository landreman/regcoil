"""`Regcoil`: assemble the current-potential regularization problem in
numpy/scipy and solve the whole lambda family from one generalized
eigendecomposition (Phase 8, ADR-021). The only Fortran call is
`regcoil._core.build_g_and_h` (Phase 7); everything else -- basis functions,
matrix assembly, the regularized solve, and the lambda scan -- is numpy/scipy.

Conventions (see docs/migration/API.md): `xn` already includes the `nfp`
factor, plasma quantities live on one field period, coil quantities span the
full torus, and the flat coil/plasma grid index used to line up with the
Fortran kernel is ``(izeta-1)*ntheta+itheta`` (itheta fastest) -- see
`_flatten_grid` below.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import scipy.linalg
import scipy.optimize

from . import _core

SYMMETRY_OPTIONS = ("stellarator_symmetric", "cos_only", "both")


def _flatten_grid(arr):
    """(..., ntheta, nzeta) -> (..., ntheta*nzeta), itheta fastest -- matches
    the Fortran kernels' flat grid index ``(izeta-1)*ntheta+itheta``."""
    return np.moveaxis(arr, -2, -1).reshape(*arr.shape[:-2], -1)


def _unflatten_grid(arr, ntheta, nzeta):
    """Inverse of `_flatten_grid`."""
    return np.moveaxis(arr.reshape(*arr.shape[:-1], nzeta, ntheta), -1, -2)


def _potential_fourier_modes(mpol, ntor, nfp):
    """Port of the legacy `regcoil_init_Fourier_modes(mpol, ntor, ...,
    include_00=.false.)`: `xm` is nonnegative, `xn` is signed (and already
    carries the `nfp` factor). The m=0 modes (n=1..ntor) come first, then
    m=1..mpol times n=-ntor..ntor -- this exact order matters because
    `single_valued_current_potential_mn` is compared element-by-element
    against legacy output in the regression tests.
    """
    xm = [0] * ntor
    xn = list(range(1, ntor + 1))
    for jm in range(1, mpol + 1):
        for jn in range(-ntor, ntor + 1):
            xm.append(jm)
            xn.append(jn)
    xm = np.array(xm, dtype=np.int64)
    xn = np.array(xn, dtype=np.int64) * nfp
    return xm, xn


class Regcoil:
    """Assembles `matrix_B`, `matrix_K`, and the RHS vectors from a `plasma`
    and `coil` surface, then caches one generalized eigendecomposition so
    every subsequent lambda is O(nbf**2) (ADR-021). Immutable after
    construction -- mutating `coil`/`plasma` afterwards does not update a
    live `Regcoil`.
    """

    def __init__(
        self,
        plasma,
        coil,
        mpol_potential,
        ntor_potential,
        net_poloidal_current=None,
        net_toroidal_current=0.0,
        symmetry="stellarator_symmetric",
    ):
        if symmetry not in SYMMETRY_OPTIONS:
            raise ValueError(f"symmetry must be one of {SYMMETRY_OPTIONS}, got {symmetry!r}")
        if coil.nfp != plasma.nfp:
            raise ValueError(f"plasma.nfp ({plasma.nfp}) != coil.nfp ({coil.nfp})")
        if net_poloidal_current is None:
            if not hasattr(plasma, "net_poloidal_current_Amperes"):
                raise ValueError(
                    "net_poloidal_current is None, but plasma has no "
                    "net_poloidal_current_Amperes attribute"
                )
            net_poloidal_current = plasma.net_poloidal_current_Amperes

        self.plasma = plasma
        self.coil = coil
        self.mpol_potential = int(mpol_potential)
        self.ntor_potential = int(ntor_potential)
        self.net_poloidal_current = float(net_poloidal_current)
        self.net_toroidal_current = float(net_toroidal_current)
        self.symmetry = symmetry
        self.nfp = nfp = plasma.nfp

        xm_potential, xn_potential = _potential_fourier_modes(self.mpol_potential, self.ntor_potential, nfp)
        mnmax_potential = xm_potential.shape[0]

        ntheta_coil, nzeta_coil = coil.ntheta, coil.nzeta
        theta_grid, zeta_grid = np.meshgrid(coil.theta, coil.zeta, indexing="ij")
        theta_flat = _flatten_grid(theta_grid)
        zeta_flat = _flatten_grid(zeta_grid)

        angle = xm_potential[:, None] * theta_flat[None, :] - xn_potential[:, None] * zeta_flat[None, :]
        sinangle = np.sin(angle)  # (mnmax_potential, ncoil_grid)
        cosangle = np.cos(angle)

        drdtheta_flat = _flatten_grid(coil.drdtheta[:, :, :nzeta_coil])  # (3, ncoil_grid)
        drdzeta_flat = _flatten_grid(coil.drdzeta[:, :, :nzeta_coil])
        norm_normal_coil_flat = _flatten_grid(coil.norm_normal)  # (ncoil_grid,)

        # Shared by f_x/f_y/f_z: xn*drdtheta + xm*drdzeta, per mode & Cartesian component.
        coef = (
            xn_potential[:, None, None] * drdtheta_flat[None, :, :]
            + xm_potential[:, None, None] * drdzeta_flat[None, :, :]
        )  # (mnmax_potential, 3, ncoil_grid)

        basis_blocks = []
        f_blocks = []
        if symmetry in ("stellarator_symmetric", "both"):
            basis_blocks.append(sinangle)
            f_blocks.append(cosangle[:, None, :] * coef)
        if symmetry in ("cos_only", "both"):
            basis_blocks.append(cosangle)
            f_blocks.append(-sinangle[:, None, :] * coef)

        if symmetry == "both":
            xm_potential = np.concatenate([xm_potential, xm_potential])
            xn_potential = np.concatenate([xn_potential, xn_potential])

        basis_all = np.concatenate(basis_blocks, axis=0)  # (nbf, ncoil_grid)
        f_all = np.concatenate(f_blocks, axis=0)  # (nbf, 3, ncoil_grid)
        nbf = basis_all.shape[0]

        d_xyz = (
            net_poloidal_current * drdtheta_flat - net_toroidal_current * drdzeta_flat
        ) / (2 * np.pi)  # (3, ncoil_grid)

        basis_functions = np.asfortranarray(basis_all.T)  # (ncoil_grid, nbf), for the Fortran kernel

        ntheta_plasma, nzeta_plasma = plasma.ntheta, plasma.nzeta
        r_plasma = np.asfortranarray(plasma.r[:, :, :nzeta_plasma])
        normal_plasma = np.asfortranarray(plasma.normal[:, :, :nzeta_plasma])
        norm_normal_plasma_flat = _flatten_grid(plasma.norm_normal)

        r_coil = np.asfortranarray(coil.r)
        normal_coil = np.asfortranarray(coil.normal)
        drdtheta_coil_all = np.asfortranarray(coil.drdtheta)
        drdzeta_coil_all = np.asfortranarray(coil.drdzeta)

        dtheta_plasma, dzeta_plasma = plasma.dtheta, plasma.dzeta
        dtheta_coil, dzeta_coil = coil.dtheta, coil.dzeta

        g, h = _core.build_g_and_h(
            r_plasma, normal_plasma, r_coil, normal_coil, drdtheta_coil_all, drdzeta_coil_all,
            basis_functions, nfp, net_poloidal_current, net_toroidal_current, dtheta_coil, dzeta_coil,
        )

        Bnormal_from_net_coil_currents = _unflatten_grid(h, ntheta_plasma, nzeta_plasma) / plasma.norm_normal
        Bnormal0 = plasma.Bnormal_from_plasma_current + Bnormal_from_net_coil_currents
        Bnormal0_flat = _flatten_grid(Bnormal0)

        RHS_B = -dtheta_plasma * dzeta_plasma * (Bnormal0_flat @ g)
        g_over_N = g / norm_normal_plasma_flat[:, None]
        matrix_B = dtheta_plasma * dzeta_plasma * (g.T @ g_over_N)

        f_over_N = f_all / norm_normal_coil_flat[None, None, :]
        matrix_K = dtheta_coil * dzeta_coil * np.einsum("mcg,ncg->mn", f_all, f_over_N)
        RHS_K = dtheta_coil * dzeta_coil * np.einsum("cg,mcg->m", d_xyz, f_over_N)

        w, V = scipy.linalg.eigh(matrix_B, matrix_K)

        # Public: the assembled problem (API.md).
        self.xm_potential = xm_potential
        self.xn_potential = xn_potential
        self.nbf = nbf
        self.basis_functions = basis_functions
        self.g = g
        self.h = h
        self.matrix_B = matrix_B
        self.matrix_K = matrix_K
        self.RHS_B = RHS_B
        self.RHS_K = RHS_K
        self.w = w
        self.V = V

        # Private: implementation detail needed to build a `Solution`.
        self._f_all = f_all
        self._d_xyz = d_xyz
        self._norm_normal_coil_flat = norm_normal_coil_flat
        self._norm_normal_plasma_flat = norm_normal_plasma_flat
        self._Bnormal0_flat = Bnormal0_flat
        self._dtheta_plasma = dtheta_plasma
        self._dzeta_plasma = dzeta_plasma
        self._dtheta_coil = dtheta_coil
        self._dzeta_coil = dzeta_coil
        self._VtRHS_B = V.T @ RHS_B
        self._VtRHS_K = V.T @ RHS_K

    def _coeffs(self, lam):
        """Mode-amplitude coefficients in the K-orthonormal eigenbasis, for a
        single (possibly infinite) lambda."""
        if np.isinf(lam):
            return self._VtRHS_K.copy()
        return (self._VtRHS_B + lam * self._VtRHS_K) / (self.w + lam)

    def solve(self, lam):
        """One regularized solve at a single lambda (`lam=np.inf` is the
        well-defined heavily-regularized limit) -> `Solution`."""
        lam = float(lam)
        solution = self.V @ self._coeffs(lam)
        return self._build_solution(lam, solution)

    def scan(self, lambdas):
        """Vectorized over an array of lambdas -- free after the
        eigendecomposition cached in `__init__`. Returns a list of
        `Solution`, one per lambda."""
        lambdas = np.atleast_1d(np.asarray(lambdas, dtype=float))
        is_inf = np.isinf(lambdas)
        finite = np.where(is_inf, 0.0, lambdas)
        coeffs = (self._VtRHS_B[:, None] + finite[None, :] * self._VtRHS_K[:, None]) / (
            self.w[:, None] + finite[None, :]
        )
        if np.any(is_inf):
            coeffs[:, is_inf] = self._VtRHS_K[:, None]
        solutions = self.V @ coeffs  # (nbf, nlambda)
        return [self._build_solution(lam, solutions[:, i]) for i, lam in enumerate(lambdas)]

    def solve_for_target(self, metric, value, xtol=1e-12, rtol=1e-12, max_iter=200):
        """Bisect (in log(lambda)) for the lambda whose `Solution.<metric>`
        (e.g. `'max_K'`, `'chi2_B'`) equals `value`. `metric` is monotonic in
        lambda; the direction is read off the `lam=0`/`lam=inf` endpoints
        rather than hard-coded, replacing the legacy staged bracket-then-
        Brent search (`regcoil_auto_regularization_solve.f90`) with a direct
        bisection on the closed-form family (ADR-021).

        Raises `ValueError` if `value` is not between the `lam=0` and
        `lam=inf` extremes (matching the legacy "target not achievable"
        exit codes -2/-3, but as an exception rather than a sentinel).
        """
        sol_lo = self.solve(0.0)
        sol_hi = self.solve(np.inf)
        f_lo = getattr(sol_lo, metric)
        f_hi = getattr(sol_hi, metric)
        achievable_lo, achievable_hi = sorted((f_lo, f_hi))
        if not (achievable_lo <= value <= achievable_hi):
            raise ValueError(
                f"target {metric}={value!r} is not achievable: the achievable range is "
                f"[{achievable_lo!r}, {achievable_hi!r}] (lambda=0 gives {f_lo!r}, "
                f"lambda=inf gives {f_hi!r})"
            )

        def residual(log_lam):
            return getattr(self.solve(np.exp(log_lam)), metric) - value

        log_lam = scipy.optimize.brentq(residual, np.log(1e-300), np.log(1e300), xtol=xtol, rtol=rtol, maxiter=max_iter)
        return self.solve(np.exp(log_lam))

    def _build_solution(self, lam, solution):
        g_sol = self.g @ solution
        Bnormal_total_flat = g_sol / self._norm_normal_plasma_flat + self._Bnormal0_flat
        Bnormal_total = _unflatten_grid(Bnormal_total_flat, self.plasma.ntheta, self.plasma.nzeta)
        max_Bnormal = float(np.max(np.abs(Bnormal_total)))
        chi2_B = float(
            self.nfp * self._dtheta_plasma * self._dzeta_plasma
            * np.sum(Bnormal_total_flat * Bnormal_total_flat * self._norm_normal_plasma_flat)
        )

        K_diff = self._d_xyz - np.einsum("mcg,m->cg", self._f_all, solution)  # (3, ncoil_grid)
        K2_times_N = np.sum(K_diff * K_diff, axis=0) / self._norm_normal_coil_flat
        chi2_K = float(self.nfp * self._dtheta_coil * self._dzeta_coil * np.sum(K2_times_N))
        max_K = float(np.sqrt(np.max(K2_times_N / self._norm_normal_coil_flat)))
        rms_K = float(np.sqrt(chi2_K / self.coil.area))

        return Solution(
            problem=self,
            lam=lam,
            solution=solution,
            chi2_B=chi2_B,
            chi2_K=chi2_K,
            max_K=max_K,
            rms_K=rms_K,
            max_Bnormal=max_Bnormal,
            Bnormal_total=Bnormal_total,
        )


@dataclass(frozen=True, eq=False)
class Solution:
    """One regularized solve, at a single lambda. `current_potential()` /
    `current_density()` are lazy (grid-sized; expanding them for every
    lambda in a scan would be wasteful).

    `eq=False`: several fields are numpy arrays, so the dataclass-generated
    `__eq__`/`__hash__` (elementwise comparison inside a bool context) would
    raise rather than compare usefully; identity equality is the sane default.
    """

    problem: Regcoil = field(repr=False)
    lam: float
    solution: np.ndarray
    chi2_B: float
    chi2_K: float
    max_K: float
    rms_K: float
    max_Bnormal: float
    Bnormal_total: np.ndarray

    @property
    def single_valued_current_potential_mn(self):
        return self.solution

    def current_potential(self):
        """(ntheta_coil, nzeta_coil) total current potential Phi (one field
        period), including the secular net-current term."""
        prob = self.problem
        coil = prob.coil
        Phi_sv = _unflatten_grid(prob.basis_functions @ self.solution, coil.ntheta, coil.nzeta)
        factor_zeta = prob.net_poloidal_current / (2 * np.pi)
        factor_theta = prob.net_toroidal_current / (2 * np.pi)
        return Phi_sv + factor_zeta * coil.zeta[None, :] + factor_theta * coil.theta[:, None]

    def current_density(self):
        """(3, ntheta_coil, nzeta_coil) Cartesian surface current density K
        (one field period)."""
        prob = self.problem
        K_diff = prob._d_xyz - np.einsum("mcg,m->cg", prob._f_all, self.solution)  # (3, ncoil_grid)
        K_flat = K_diff / prob._norm_normal_coil_flat[None, :]
        return _unflatten_grid(K_flat, prob.coil.ntheta, prob.coil.nzeta)
