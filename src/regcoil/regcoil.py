"""`Regcoil`: assemble the current-potential regularization problem in
numpy/scipy and solve the whole lambda family from one generalized
eigendecomposition. The only Fortran call is
`regcoil._core.build_g_and_h`. Everything else -- basis functions,
matrix assembly, the regularized solve, and the lambda scan -- is numpy/scipy.

Conventions (see docs/migration/API.md): `xn` already includes the `nfp`
factor, plasma quantities live on one field period, coil quantities span the
full torus, and the flat coil/plasma grid index used to line up with the
Fortran kernel is ``(izeta-1)*ntheta+itheta`` (itheta fastest) -- see
`_flatten_grid` below.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass, field
import logging
from time import perf_counter
from typing import cast

import numpy as np
import scipy.linalg
import scipy.linalg.blas
import scipy.optimize

from . import _core
from .magnetic_field import (
    DEFAULT_CHUNK_SIZE as _MAGNETIC_FIELD_CHUNK_SIZE,
    MagneticFieldEvaluator,
)

logger = logging.getLogger(__name__)

_F_ALL_MODE_BLOCK = 64
_MATRIX_K_GRID_BLOCK = 2048


def _flatten_grid(arr):
    """(..., ntheta, nzeta) -> (..., ntheta*nzeta), itheta fastest -- matches
    the Fortran kernels' flat grid index ``(izeta-1)*ntheta+itheta``."""
    return np.moveaxis(arr, -2, -1).reshape(*arr.shape[:-2], -1)


def _unflatten_grid(arr, ntheta, nzeta):
    """Inverse of `_flatten_grid`."""
    return np.moveaxis(arr.reshape(*arr.shape[:-1], nzeta, ntheta), -1, -2)


def _build_matrix_K(f_all, norm_normal_coil_flat, scale, grid_block=_MATRIX_K_GRID_BLOCK):
    """Assemble the symmetric K Gram matrix without a full weighted copy of f_all."""
    if grid_block < 1:
        raise ValueError(f"grid_block must be positive, got {grid_block}")

    nbf, _, ncoil_grid = f_all.shape
    sqrtN_coil_flat = np.sqrt(norm_normal_coil_flat)
    K_upper = np.zeros((nbf, nbf), dtype=np.float64, order="F")

    first_block = True
    for start in range(0, ncoil_grid, grid_block):
        stop = min(start + grid_block, ncoil_grid)
        weighted_block = (
            f_all[:, :, start:stop]
            / sqrtN_coil_flat[None, None, start:stop]
        )
        # weighted_block is C-contiguous, so this transpose is Fortran-contiguous
        # without another copy. DSYRK then accumulates A.T @ A into one triangle.
        A_block = np.asfortranarray(weighted_block.reshape(nbf, -1).T)
        K_upper = scipy.linalg.blas.dsyrk(
            alpha=scale,
            a=A_block,
            trans=1,
            beta=0.0 if first_block else 1.0,
            c=K_upper,
            lower=0,
            overwrite_c=1,
        )
        first_block = False

    return np.triu(K_upper) + np.triu(K_upper, 1).T


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
    every subsequent lambda is O(nbf**2). Immutable after
    construction -- mutating `coil`/`plasma` afterwards does not update a
    live `Regcoil`.

    `__init__` is split into a cheap assembly (surfaces, potential modes,
    `basis_functions`, `f_all`, `d_xyz` -- no Fortran, no big BLAS) and an
    expensive one (`_build_operators`: `build_g_and_h`, `matrix_B/K`,
    `eigh`), so that `regcoil.load()` can reconstruct a `Regcoil`
    from disk running only the cheap part -- its stored Solutions then plot
    with no kernel call. The operators are rebuilt lazily, via
    `_ensure_operators`, only if a *new* lambda is solved on such an object.
    """

    def __init__(
        self,
        plasma,
        coil,
        mpol_potential,
        ntor_potential,
        net_poloidal_current=None,
        net_toroidal_current=0.0,
        stellarator_symmetric=True,
    ):
        self._init_cheap(
            plasma, coil, mpol_potential, ntor_potential,
            net_poloidal_current, net_toroidal_current, stellarator_symmetric,
        )
        self._build_operators()

    @classmethod
    def _from_loaded(
        cls, plasma, coil, mpol_potential, ntor_potential,
        net_poloidal_current, net_toroidal_current, stellarator_symmetric,
        Bnormal_from_net_coil_currents,
    ):
        """Used by `regcoil.load()`: cheap assembly only. The operators
        (`g`, `h`, `matrix_B/K`, `w`, `V`) are left unset and rebuilt lazily
        (see `_ensure_operators`) only if `solve()`/`scan()` is later called.
        `Bnormal_from_net_coil_currents` is the one operator-derived
        quantity stored on disk (its recompute needs the Fortran `h` term),
        so it is set here directly rather than recomputed.
        """
        obj = cls.__new__(cls)
        obj._init_cheap(
            plasma, coil, mpol_potential, ntor_potential,
            net_poloidal_current, net_toroidal_current, stellarator_symmetric,
        )
        obj.Bnormal_from_net_coil_currents = Bnormal_from_net_coil_currents
        return obj

    def _init_cheap(
        self, plasma, coil, mpol_potential, ntor_potential,
        net_poloidal_current, net_toroidal_current, stellarator_symmetric,
    ):
        if not isinstance(stellarator_symmetric, (bool, np.bool_)):
            raise TypeError(
                "stellarator_symmetric must be a bool: "
                "True for stellarator-symmetric basis only, "
                "False for both sine and cosine basis families"
            )
        if coil.nfp != plasma.nfp:
            raise ValueError(f"plasma.nfp ({plasma.nfp}) != coil.nfp ({coil.nfp})")
        if net_poloidal_current is None:
            if not hasattr(plasma, "net_poloidal_current"):
                raise ValueError(
                    "net_poloidal_current is None, but plasma has no "
                    "net_poloidal_current attribute"
                )
            net_poloidal_current = plasma.net_poloidal_current

        self.plasma = plasma
        self.coil = coil
        self.mpol_potential = int(mpol_potential)
        self.ntor_potential = int(ntor_potential)
        self.net_poloidal_current = float(net_poloidal_current)
        self.net_toroidal_current = float(net_toroidal_current)
        self.stellarator_symmetric = bool(stellarator_symmetric)
        self.nfp = nfp = plasma.nfp

        xm_potential, xn_potential = _potential_fourier_modes(self.mpol_potential, self.ntor_potential, nfp)

        nzeta_coil = coil.nzeta
        theta_grid, zeta_grid = np.meshgrid(coil.theta, coil.zeta, indexing="ij")
        theta_flat = _flatten_grid(theta_grid)
        zeta_flat = _flatten_grid(zeta_grid)

        angle = xm_potential[:, None] * theta_flat[None, :] - xn_potential[:, None] * zeta_flat[None, :]
        sinangle = np.sin(angle)  # (mnmax_potential, ncoil_grid)
        np.cos(angle, out=angle)
        cosangle = angle

        drdtheta_flat = _flatten_grid(coil.drdtheta[:, :, :nzeta_coil])  # (3, ncoil_grid)
        drdzeta_flat = _flatten_grid(coil.drdzeta[:, :, :nzeta_coil])
        norm_normal_coil_flat = _flatten_grid(coil.norm_normal)  # (ncoil_grid,)

        if self.stellarator_symmetric:
            # Building the full coef and then concatenating one-element lists
            # used three simultaneous nmode*3*ngrid arrays. Fill the final
            # allocation in mode blocks instead.
            basis_all = sinangle
            nbf = basis_all.shape[0]
            f_all = np.empty((nbf, 3, theta_flat.size), dtype=np.float64)
            for start in range(0, nbf, _F_ALL_MODE_BLOCK):
                stop = min(start + _F_ALL_MODE_BLOCK, nbf)
                sl = slice(start, stop)
                coef_block = (
                    xn_potential[sl, None, None] * drdtheta_flat[None, :, :]
                    + xm_potential[sl, None, None] * drdzeta_flat[None, :, :]
                )
                np.multiply(cosangle[sl, None, :], coef_block, out=f_all[sl])
                del coef_block
        else:
            # Shared by f_x/f_y/f_z: xn*drdtheta + xm*drdzeta, per mode
            # and Cartesian component.
            coef = (
                xn_potential[:, None, None] * drdtheta_flat[None, :, :]
                + xm_potential[:, None, None] * drdzeta_flat[None, :, :]
            )
            basis_all = np.concatenate([sinangle, cosangle], axis=0)
            f_all = np.concatenate(
                [cosangle[:, None, :] * coef, -sinangle[:, None, :] * coef],
                axis=0,
            )
            xm_potential = np.concatenate([xm_potential, xm_potential])
            xn_potential = np.concatenate([xn_potential, xn_potential])
            nbf = basis_all.shape[0]

        d_xyz = (
            net_poloidal_current * drdtheta_flat - net_toroidal_current * drdzeta_flat
        ) / (2 * np.pi)  # (3, ncoil_grid)

        basis_functions = np.asfortranarray(basis_all.T)  # (ncoil_grid, nbf), for the Fortran kernel

        norm_normal_plasma_flat = _flatten_grid(plasma.norm_normal)
        dtheta_plasma, dzeta_plasma = plasma.dtheta, plasma.dzeta
        dtheta_coil, dzeta_coil = coil.dtheta, coil.dzeta

        # Public: the assembled problem (API.md).
        self.xm_potential = xm_potential
        self.xn_potential = xn_potential
        self.nbf = nbf
        self.basis_functions = basis_functions

        # Private: implementation detail needed to build a `Solution`.
        self._f_all = f_all
        self._d_xyz = d_xyz
        self._norm_normal_coil_flat = norm_normal_coil_flat
        self._norm_normal_plasma_flat = norm_normal_plasma_flat
        self._dtheta_plasma = dtheta_plasma
        self._dzeta_plasma = dzeta_plasma
        self._dtheta_coil = dtheta_coil
        self._dzeta_coil = dzeta_coil

        # Operators (expensive assembly, `_build_operators`) -- unset until built.
        self.g = self.h = None
        self.matrix_B = self.matrix_K = None
        self.RHS_B = self.RHS_K = None
        self.w = self.V = None
        self.Bnormal_from_net_coil_currents = None
        self._Bnormal0_flat = None
        self._VtRHS_B = self._VtRHS_K = None

    def _ensure_operators(self):
        """Lazily rebuild the operators (Fortran kernel, big matrices,
        eigendecomposition) if this `Regcoil` came from `regcoil.load()`
        without them."""
        if self.g is None:
            logger.info("Rebuilding operators for a loaded Regcoil (new lambda solve requested)")
            self._build_operators()

    def _build_operators(self):
        plasma, coil = self.plasma, self.coil
        nfp = self.nfp
        net_poloidal_current = self.net_poloidal_current
        net_toroidal_current = self.net_toroidal_current
        basis_functions = self.basis_functions
        f_all = self._f_all
        d_xyz = self._d_xyz
        norm_normal_coil_flat = self._norm_normal_coil_flat
        norm_normal_plasma_flat = self._norm_normal_plasma_flat
        dtheta_plasma, dzeta_plasma = self._dtheta_plasma, self._dzeta_plasma
        dtheta_coil, dzeta_coil = self._dtheta_coil, self._dzeta_coil
        nbf = self.nbf

        ntheta_plasma, nzeta_plasma = plasma.ntheta, plasma.nzeta
        ntheta_coil, nzeta_coil = coil.ntheta, coil.nzeta
        r_plasma = np.asfortranarray(plasma.r[:, :, :nzeta_plasma])
        normal_plasma = np.asfortranarray(plasma.normal[:, :, :nzeta_plasma])

        r_coil = np.asfortranarray(coil.r)
        normal_coil = np.asfortranarray(coil.normal)
        drdtheta_coil_all = np.asfortranarray(coil.drdtheta)
        drdzeta_coil_all = np.asfortranarray(coil.drdzeta)

        logger.info(
            "Starting build_g_and_h for plasma %dx%d, coil %dx%d, nbf=%d",
            ntheta_plasma,
            nzeta_plasma,
            ntheta_coil,
            nzeta_coil,
            nbf,
        )
        kernel_start = perf_counter()
        g, h = _core.build_g_and_h(  # type: ignore[attr-defined]
            r_plasma, normal_plasma, r_coil, normal_coil, drdtheta_coil_all, drdzeta_coil_all,
            basis_functions, nfp, net_poloidal_current, net_toroidal_current, dtheta_coil, dzeta_coil,
        )
        logger.info(
            "Finished build_g_and_h in %.3f s",
            perf_counter() - kernel_start,
        )

        logger.info("Starting Bnormal0_flat assembly")
        assembly_start = perf_counter()
        Bnormal_from_net_coil_currents = _unflatten_grid(h, ntheta_plasma, nzeta_plasma) / plasma.norm_normal
        Bnormal0 = plasma.Bnormal_from_plasma_current + Bnormal_from_net_coil_currents
        Bnormal0_flat = _flatten_grid(Bnormal0)
        logger.info(
            "Finished Bnormal0_flat assembly in %.3f s",
            perf_counter() - assembly_start,
        )

        logger.info("Starting RHS_B assembly")
        start_time = perf_counter()
        RHS_B = -dtheta_plasma * dzeta_plasma * (Bnormal0_flat @ g)
        logger.info("Finished RHS_B assembly in %.3f s", perf_counter() - start_time)
        g_over_N = g / norm_normal_plasma_flat[:, None]
        logger.info("Starting matrix_B assembly")
        start_time = perf_counter()
        matrix_B = dtheta_plasma * dzeta_plasma * (g.T @ g_over_N)
        logger.info("Finished matrix_B assembly in %.3f s", perf_counter() - start_time)

        logger.info("Starting matrix_K assembly")
        start_time = perf_counter()
        matrix_K = _build_matrix_K(
            f_all,
            norm_normal_coil_flat,
            dtheta_coil * dzeta_coil,
        )
        logger.info("Finished matrix_K assembly in %.3f s", perf_counter() - start_time)
        logger.info("Starting RHS_K assembly")
        start_time = perf_counter()
        weighted_d_xyz = d_xyz / norm_normal_coil_flat[None, :]
        RHS_K = dtheta_coil * dzeta_coil * np.einsum(
            "cg,mcg->m", weighted_d_xyz, f_all, optimize=True,
        )
        logger.info("Finished RHS_K assembly in %.3f s", perf_counter() - start_time)

        logger.info("Starting generalized eigensolve for %d basis functions", nbf)
        eigensolve_start = perf_counter()
        w, V = scipy.linalg.eigh(matrix_B, matrix_K)
        logger.info(
            "Finished generalized eigensolve in %.3f s",
            perf_counter() - eigensolve_start,
        )

        self.g = g
        self.h = h
        self.matrix_B = matrix_B
        self.matrix_K = matrix_K
        self.RHS_B = RHS_B
        self.RHS_K = RHS_K
        self.w = w
        self.V = V
        self.Bnormal_from_net_coil_currents = Bnormal_from_net_coil_currents
        self._Bnormal0_flat = Bnormal0_flat
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
        self._ensure_operators()
        lam = float(lam)
        solution = self.V @ self._coeffs(lam)
        return self._build_solution(lam, solution)

    def scan(self, lambdas):
        """Vectorized over an array of lambdas -- free after the
        eigendecomposition cached in `__init__`. Returns a `SolutionScan`
        (a `Sequence[Solution]` with columnar `.lam`/`.f_B`/... accessors)."""
        self._ensure_operators()
        logger.info("Starting lambda scan")
        start_time = perf_counter()
        lambdas = np.atleast_1d(np.asarray(lambdas, dtype=float))
        is_inf = np.isinf(lambdas)
        finite = np.where(is_inf, 0.0, lambdas)
        coeffs = (self._VtRHS_B[:, None] + finite[None, :] * self._VtRHS_K[:, None]) / (
            self.w[:, None] + finite[None, :]
        )
        if np.any(is_inf):
            coeffs[:, is_inf] = self._VtRHS_K[:, None]
        solutions = self.V @ coeffs  # (nbf, nlambda)
        data = [self._build_solution(lam, solutions[:, i]) for i, lam in enumerate(lambdas)]
        logger.info(
            "Finished lambda scan in %.3f s",
            perf_counter() - start_time,
        )
        return SolutionScan(data)

    def solve_for_target(self, metric, value, xtol=1e-12, rtol=1e-12, max_iter=200):
        """Bisect (in log(lambda)) for the lambda whose `Solution.<metric>`
        (e.g. `'max_K'`, `'f_B'`) equals `value`. `metric` is monotonic in
        lambda; the direction is read off the `lam=0`/`lam=inf` endpoints
        rather than hard-coded, replacing the legacy staged bracket-then-
        Brent search (`regcoil_auto_regularization_solve.f90`) with a direct
        bisection on the closed-form family.

        Raises `ValueError` if `value` is not between the `lam=0` and
        `lam=inf` extremes (matching the legacy "target not achievable"
        exit codes -2/-3, but as an exception rather than a sentinel).
        """
        self._ensure_operators()
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

        logger.info("Starting lambda search")
        start_time = perf_counter()
        log_lam = cast(
            float,
            scipy.optimize.brentq(  # pyright: ignore[reportArgumentType,reportCallIssue]
                residual,
                np.log(1e-300),
                np.log(1e300),
                xtol=xtol,
                rtol=rtol,  # pyright: ignore[reportArgumentType]
                maxiter=max_iter,
            ),
        )
        logger.info(
            "Finished lambda search in %.3f s",
            perf_counter() - start_time,
        )
        return self.solve(np.exp(log_lam))

    def _build_solution(self, lam, solution):
        logger.info("Starting solution build for lambda=%g", lam)
        start_time = perf_counter()
        g_sol = self.g @ solution
        Bnormal_total_flat = g_sol / self._norm_normal_plasma_flat + self._Bnormal0_flat
        Bnormal_total = _unflatten_grid(Bnormal_total_flat, self.plasma.ntheta, self.plasma.nzeta)
        max_Bnormal = float(np.max(np.abs(Bnormal_total)))
        f_B = float(
            self.nfp * self._dtheta_plasma * self._dzeta_plasma
            * np.sum(Bnormal_total_flat * Bnormal_total_flat * self._norm_normal_plasma_flat)
        )
        abs_Bnormal_over_B = np.abs(Bnormal_total / self.plasma.modB)
        max_Bnormal_over_B = float(np.max(abs_Bnormal_over_B))
        avg_Bnormal_over_B = float(
            np.sum(abs_Bnormal_over_B * self.plasma.norm_normal)
            / np.sum(self.plasma.norm_normal)
        )

        K_diff = self._d_xyz - np.einsum("mcg,m->cg", self._f_all, solution)  # (3, ncoil_grid)
        K2_times_N = np.sum(K_diff * K_diff, axis=0) / self._norm_normal_coil_flat
        f_K = float(self.nfp * self._dtheta_coil * self._dzeta_coil * np.sum(K2_times_N))
        max_K = float(np.sqrt(np.max(K2_times_N / self._norm_normal_coil_flat)))
        rms_K = float(np.sqrt(f_K / self.coil.area))
        logger.info(
            "Finished solution build for lambda=%g in %.3f s",
            lam,
            perf_counter() - start_time,
        )

        return Solution(
            problem=self,
            lam=lam,
            solution=solution,
            f_B=f_B,
            f_K=f_K,
            max_K=max_K,
            rms_K=rms_K,
            max_Bnormal=max_Bnormal,
            max_Bnormal_over_B=max_Bnormal_over_B,
            avg_Bnormal_over_B=avg_Bnormal_over_B,
            Bnormal_total=Bnormal_total,
        )

    def save(self, path):
        """Save this problem (plus its plasma and coil) to `path`."""
        from . import _serialize
        _serialize.save(path, problem=self)

    @classmethod
    def load(cls, path):
        """Load a `Regcoil` (plus its plasma and coil) from `path`. The
        expensive operators are rebuilt lazily only if a new lambda is
        solved."""
        from . import _serialize
        return _serialize.load_problem(path)


@dataclass(frozen=True, eq=False)
class Solution:
    """One regularized solve, at a single lambda. `current_potential()` /
    `current_density()` are lazy (grid-sized; expanding them for every
    lambda in a scan would be wasteful) for a freshly-computed `Solution`,
    but a `Solution` reconstructed by `regcoil.load()` carries the
    precomputed grids (`_current_potential`/`_current_density`, stored on
    disk to avoid materializing the ~0.5 GB `f_all` intermediate on load),
    so those accessors are O(1) for a loaded run.

    `eq=False`: several fields are numpy arrays, so the dataclass-generated
    `__eq__`/`__hash__` (elementwise comparison inside a bool context) would
    raise rather than compare usefully; identity equality is the sane default.
    """

    problem: Regcoil = field(repr=False)
    lam: float
    solution: np.ndarray
    f_B: float
    f_K: float
    max_K: float
    rms_K: float
    max_Bnormal: float
    max_Bnormal_over_B: float
    avg_Bnormal_over_B: float
    Bnormal_total: np.ndarray
    _current_potential: np.ndarray | None = field(default=None, repr=False)
    _current_density: np.ndarray | None = field(default=None, repr=False)
    _magnetic_field_evaluators: dict[
        tuple[int, int], MagneticFieldEvaluator
    ] = field(
        default_factory=dict, init=False, repr=False, compare=False
    )

    @property
    def single_valued_current_potential_mn(self):
        return self.solution

    def current_potential(self):
        """(ntheta_coil, nzeta_coil) total current potential Phi (one field
        period), including the secular net-current term."""
        if self._current_potential is not None:
            return self._current_potential
        prob = self.problem
        coil = prob.coil
        Phi_sv = _unflatten_grid(prob.basis_functions @ self.solution, coil.ntheta, coil.nzeta)
        factor_zeta = prob.net_poloidal_current / (2 * np.pi)
        factor_theta = prob.net_toroidal_current / (2 * np.pi)
        return Phi_sv + factor_zeta * coil.zeta[None, :] + factor_theta * coil.theta[:, None]

    def current_density(self):
        """(3, ntheta_coil, nzeta_coil) Cartesian surface current density K
        (one field period)."""
        if self._current_density is not None:
            return self._current_density
        prob = self.problem
        K_diff = prob._d_xyz - np.einsum("mcg,m->cg", prob._f_all, self.solution)  # (3, ncoil_grid)
        K_flat = K_diff / prob._norm_normal_coil_flat[None, :]
        return _unflatten_grid(K_flat, prob.coil.ntheta, prob.coil.nzeta)

    def prepare_magnetic_field(self, *, theta_stride=1, zeta_stride=1):
        """Return a cached, prepared Biot--Savart evaluator.

        ``theta_stride`` subsamples the one-period poloidal grid and
        ``zeta_stride`` subsamples the assembled full-torus toroidal grid.
        Each must divide its corresponding periodic grid size. The retained
        nodes are reweighted as a uniform trapezoidal quadrature. Values above
        one trade quadrature accuracy for field-evaluation speed and should be
        selected through a convergence study for the intended diagnostics.
        The returned evaluator's ``quadrature_shape`` is ordered as
        ``(full_torus_zeta, theta)``.
        """
        theta_stride, zeta_stride = MagneticFieldEvaluator._validated_strides(
            self, theta_stride, zeta_stride
        )
        key = (theta_stride, zeta_stride)
        evaluator = self._magnetic_field_evaluators.get(key)
        if evaluator is None:
            evaluator = MagneticFieldEvaluator.from_solution(
                self,
                theta_stride=theta_stride,
                zeta_stride=zeta_stride,
            )
            self._magnetic_field_evaluators[key] = evaluator
        return evaluator

    def magnetic_field(
        self,
        points,
        *,
        chunk_size=_MAGNETIC_FIELD_CHUNK_SIZE,
        theta_stride=1,
        zeta_stride=1,
    ):
        """Magnetic field produced by the coil-surface current.

        Parameters
        ----------
        points : array_like, shape (..., 3)
          Cartesian evaluation points in meters.
        chunk_size : int
          Number of evaluation points processed simultaneously.
        theta_stride, zeta_stride : int
          Uniform quadrature subsampling factors. See
          :meth:`prepare_magnetic_field`.

        Returns
        -------
        B : ndarray, shape (..., 3)
          Cartesian magnetic field in tesla.

        Notes
        -----
        This computes the field of the continuous REGCOIL surface-current
        solution. It does not include the plasma-current field and does not
        compute the field of discrete coils produced by ``regcoil.cut``.

        Evaluation on or extremely close to the coil surface is singular and
        is not supported by this direct quadrature.
        """
        evaluator = self.prepare_magnetic_field(
            theta_stride=theta_stride, zeta_stride=zeta_stride
        )
        return evaluator(points, chunk_size=chunk_size)

    def save(self, path):
        """Save this solution (plus its problem, plasma, and coil) to `path`."""
        from . import _serialize
        _serialize.save(path, solutions=[self])

    @classmethod
    def load(cls, path):
        """Load a single-lambda `Solution` (plus its problem, plasma, and
        coil) from `path`. Raises `ValueError` if `path` holds more than one
        lambda -- use `regcoil.load(path).solutions` for a scan."""
        from . import _serialize
        return _serialize.load_solution(path)

    def plot_current_potential(self, kind="single_valued", ax=None):
        """Convenience delegate for `regcoil.plot.current_potential`."""
        from . import plot
        return plot.current_potential(self, kind=kind, ax=ax)

    def plot_current_density(self, ax=None):
        """Convenience delegate for `regcoil.plot.current_density`."""
        from . import plot
        return plot.current_density(self, ax=ax)

    def plot_bnormal(self, component="total", ax=None):
        """Convenience delegate for `regcoil.plot.bnormal`."""
        from . import plot
        return plot.bnormal(self, component=component, ax=ax)

    def cut(self, coils_per_half_period, thickness=None, theta_shift=0):
        """Convenience delegate for `regcoil.cut.cut`."""
        from . import cut as cut_module
        return cut_module.cut(
            self, coils_per_half_period, thickness=thickness, theta_shift=theta_shift
        )


class SolutionScan(Sequence):
    """A lambda scan: a `Sequence[Solution]` that iterates/indexes as
    ordinary `Solution` objects, plus columnar `.lam`/`.f_B`/`.f_K`/`.max_K`/
    `.rms_K`/`.max_Bnormal`/`.max_Bnormal_over_B`/`.avg_Bnormal_over_B` array
    accessors for direct Pareto/scan plotting. Returned by `Regcoil.scan()`
    and `regcoil.load(...).solutions`.
    """

    def __init__(self, solutions):
        self._solutions = list(solutions)

    def __len__(self):
        return len(self._solutions)

    def __getitem__(self, index):
        if isinstance(index, slice):
            return SolutionScan(self._solutions[index])
        return self._solutions[index]

    def __repr__(self):
        return f"SolutionScan({len(self._solutions)} solutions)"

    def _column(self, name):
        return np.array([getattr(sol, name) for sol in self._solutions])

    @property
    def lam(self):
        return self._column("lam")

    @property
    def f_B(self):
        return self._column("f_B")

    @property
    def f_K(self):
        return self._column("f_K")

    @property
    def max_K(self):
        return self._column("max_K")

    @property
    def rms_K(self):
        return self._column("rms_K")

    @property
    def max_Bnormal(self):
        return self._column("max_Bnormal")

    @property
    def max_Bnormal_over_B(self):
        return self._column("max_Bnormal_over_B")

    @property
    def avg_Bnormal_over_B(self):
        return self._column("avg_Bnormal_over_B")

    def plot_pareto(self, x="f_K", y="f_B", ax=None):
        """Convenience delegate for `regcoil.plot.pareto`."""
        from . import plot
        return plot.pareto(self, x=x, y=y, ax=ax)

    def plot_lambda_scan(self, ax=None):
        """Convenience delegate for `regcoil.plot.lambda_scan`."""
        from . import plot
        return plot.lambda_scan(self, ax=ax)

    def plot_current_potential_scan(self, kind="total", nmax=16, fig=None):
        """Convenience delegate for `regcoil.plot.current_potential_scan`."""
        from . import plot
        return plot.current_potential_scan(self, kind=kind, nmax=nmax, fig=fig)

    def plot_bnormal_scan(self, nmax=16, fig=None):
        """Convenience delegate for `regcoil.plot.bnormal_scan`."""
        from . import plot
        return plot.bnormal_scan(self, nmax=nmax, fig=fig)
