"""Toroidal-angle field-line tracing and derived diagnostics.

The routines in this module are independent of plotting and accept any
callable magnetic field returning Cartesian ``(B_x, B_y, B_z)``. Use
:meth:`ToroidalTracer.from_solution` for REGCOIL's prepared surface-current
field.
"""

from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
import math
import operator
from time import perf_counter
from typing import Callable, Sequence

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import least_squares


Array = np.ndarray
MagneticField = Callable[[Array], Array]


@dataclass(frozen=True)
class ToroidalTracingOptions:
    """ODE and conditioning options for toroidal-angle tracing."""

    bphi_absolute_floor: float = 1.0e-10
    bphi_relative_floor: float = 1.0e-3
    maximum_rhs_norm: float = 100.0
    minimum_field_magnitude: float = 1.0e-14
    method: str = "DOP853"
    rtol: float = 1.0e-8
    atol: float = 1.0e-10

    def __post_init__(self):
        positive = (
            "bphi_absolute_floor",
            "bphi_relative_floor",
            "maximum_rhs_norm",
            "minimum_field_magnitude",
            "rtol",
            "atol",
        )
        for name in positive:
            value = getattr(self, name)
            if not np.isfinite(value) or value <= 0:
                raise ValueError(f"{name} must be finite and positive")


@dataclass(frozen=True)
class ToroidalDiagnostic:
    """Conditioning information at one cylindrical field-line position."""

    valid: bool
    reason: str
    phi: float
    R: float
    Z: float
    B_R: float = math.nan
    B_phi: float = math.nan
    B_Z: float = math.nan
    B_magnitude: float = math.nan
    bphi_fraction: float = math.nan
    rhs_norm: float = math.nan
    rhs: Array | None = None


class ToroidalTracingError(RuntimeError):
    """A field line cannot safely use toroidal angle as its parameter."""

    def __init__(self, diagnostic: ToroidalDiagnostic):
        self.diagnostic = diagnostic
        super().__init__(diagnostic.reason)


@dataclass(frozen=True)
class PoincareTrace:
    """One successful or expected-failure toroidal field-line trace."""

    initial_rz: Array
    success: bool
    phi: Array
    rz: Array
    nfev: int
    message: str
    failure_diagnostic: ToroidalDiagnostic | None
    minimum_abs_B_phi: float
    minimum_bphi_fraction: float
    maximum_rhs_norm: float
    phase_fractions: Array
    phase_sample_indices: tuple[Array, ...]
    direction: int
    nfp: int
    periods: int

    def section(self, phase_index: int, *, omit_initial=False):
        """Return ``(phi, rz)`` samples for one requested physical section."""
        phase_index = operator.index(phase_index)
        indices = self.phase_sample_indices[phase_index]
        indices = indices[indices < len(self.phi)]
        if (
            omit_initial
            and len(indices)
            and indices[0] == 0
            and np.isclose(self.phase_fractions[phase_index], 0.0)
        ):
            indices = indices[1:]
        return self.phi[indices], self.rz[indices]

    def same_section_returns(self):
        """Return samples at repeated returns to the physical phi=0 section."""
        matches = np.flatnonzero(np.isclose(self.phase_fractions, 0.0))
        if len(matches) != 1:
            raise ValueError("phase_fractions must contain phi=0 exactly once")
        return self.section(int(matches[0]))


@dataclass(frozen=True)
class PoincareBatch:
    """Ordered results and shared metadata for an independent trace batch."""

    traces: tuple[PoincareTrace, ...]
    phase_fractions: Array
    direction: int
    nfp: int
    periods: int
    wall_time: float

    @property
    def successful(self):
        return tuple(trace for trace in self.traces if trace.success)

    @property
    def failed(self):
        return tuple(trace for trace in self.traces if not trace.success)


@dataclass(frozen=True)
class MagneticAxisResult:
    """Accepted elliptic fixed point of a one-field-period return map."""

    rz: Array
    initial_guess: Array
    return_residual: Array
    return_error: float
    return_map_jacobian: Array
    eigenvalues: Array
    determinant: float
    discriminant: float
    return_map_evaluations: int
    optimizer_message: str
    wall_time: float


@dataclass(frozen=True)
class IotaProfile:
    """Geometric rotational-transform estimates ordered by radius."""

    radius: Array
    initial_rz: Array
    period_counts: Array
    iota: Array
    axis_rz: Array


def cylindrical_field(phi, rz, magnetic_field: MagneticField):
    """Evaluate ``(B_R, B_phi, B_Z)`` at cylindrical ``(R, phi, Z)``."""
    R, Z = np.asarray(rz, dtype=float)
    cosine = np.cos(phi)
    sine = np.sin(phi)
    B_x, B_y, B_Z = magnetic_field(
        np.array([R * cosine, R * sine, Z])
    )
    B_R = cosine * B_x + sine * B_y
    B_phi = -sine * B_x + cosine * B_y
    return float(B_R), float(B_phi), float(B_Z)


def _invalid(reason, phi, R, Z, **values):
    diagnostic = ToroidalDiagnostic(
        valid=False,
        reason=reason,
        phi=float(phi),
        R=float(R),
        Z=float(Z),
        **values,
    )
    raise ToroidalTracingError(diagnostic)


def _evaluate_toroidal_rhs(phi, rz, magnetic_field, options):
    """Hot scalar RHS path; diagnostics are allocated only on failure."""
    R = float(rz[0])
    Z = float(rz[1])
    if not np.isfinite(R) or not np.isfinite(Z) or R <= 0:
        _invalid(f"invalid cylindrical position (R, Z)=({R}, {Z})", phi, R, Z)

    B_R, B_phi, B_Z = cylindrical_field(phi, (R, Z), magnetic_field)
    if not (np.isfinite(B_R) and np.isfinite(B_phi) and np.isfinite(B_Z)):
        _invalid(
            "magnetic field contains a non-finite component",
            phi,
            R,
            Z,
            B_R=B_R,
            B_phi=B_phi,
            B_Z=B_Z,
        )

    magnitude = float(np.sqrt(B_R * B_R + B_phi * B_phi + B_Z * B_Z))
    if magnitude < options.minimum_field_magnitude:
        _invalid(
            f"magnetic-field magnitude is too small: |B|={magnitude:.3e} T",
            phi,
            R,
            Z,
            B_R=B_R,
            B_phi=B_phi,
            B_Z=B_Z,
            B_magnitude=magnitude,
        )

    absolute_B_phi = abs(B_phi)
    bphi_fraction = absolute_B_phi / magnitude
    if absolute_B_phi < options.bphi_absolute_floor:
        _invalid(
            f"|B_phi|={absolute_B_phi:.3e} T is below the absolute floor "
            f"{options.bphi_absolute_floor:.3e} T",
            phi,
            R,
            Z,
            B_R=B_R,
            B_phi=B_phi,
            B_Z=B_Z,
            B_magnitude=magnitude,
            bphi_fraction=bphi_fraction,
        )

    rhs_R = R * B_R / B_phi
    rhs_Z = R * B_Z / B_phi
    rhs_norm = float(math.hypot(rhs_R, rhs_Z))
    if not (np.isfinite(rhs_R) and np.isfinite(rhs_Z)):
        _invalid(
            "the cylindrical field-line RHS is non-finite",
            phi,
            R,
            Z,
            B_R=B_R,
            B_phi=B_phi,
            B_Z=B_Z,
            B_magnitude=magnitude,
            bphi_fraction=bphi_fraction,
            rhs_norm=rhs_norm,
            rhs=np.array([rhs_R, rhs_Z]),
        )
    if bphi_fraction < options.bphi_relative_floor:
        _invalid(
            f"|B_phi|/|B|={bphi_fraction:.3e} is below the relative floor "
            f"{options.bphi_relative_floor:.3e}",
            phi,
            R,
            Z,
            B_R=B_R,
            B_phi=B_phi,
            B_Z=B_Z,
            B_magnitude=magnitude,
            bphi_fraction=bphi_fraction,
            rhs_norm=rhs_norm,
            rhs=np.array([rhs_R, rhs_Z]),
        )
    if rhs_norm > options.maximum_rhs_norm:
        _invalid(
            f"|d(R,Z)/dphi|={rhs_norm:.3e} m/rad exceeds "
            f"{options.maximum_rhs_norm:.3e} m/rad",
            phi,
            R,
            Z,
            B_R=B_R,
            B_phi=B_phi,
            B_Z=B_Z,
            B_magnitude=magnitude,
            bphi_fraction=bphi_fraction,
            rhs_norm=rhs_norm,
            rhs=np.array([rhs_R, rhs_Z]),
        )
    return (
        np.array([rhs_R, rhs_Z]),
        absolute_B_phi,
        bphi_fraction,
        rhs_norm,
    )


class _TraceMonitor:
    def __init__(self, tracer, options):
        self.tracer = tracer
        self.options = options
        self.minimum_abs_B_phi = np.inf
        self.minimum_bphi_fraction = np.inf
        self.maximum_rhs_norm = 0.0
        self.last_diagnostic = None

    def __call__(self, phi, rz):
        try:
            rhs, absolute, fraction, norm = _evaluate_toroidal_rhs(
                phi, rz, self.tracer.magnetic_field, self.options
            )
        except ToroidalTracingError as error:
            self.last_diagnostic = error.diagnostic
            diagnostic = error.diagnostic
            if np.isfinite(diagnostic.B_phi):
                self.minimum_abs_B_phi = min(
                    self.minimum_abs_B_phi, abs(diagnostic.B_phi)
                )
            if np.isfinite(diagnostic.bphi_fraction):
                self.minimum_bphi_fraction = min(
                    self.minimum_bphi_fraction, diagnostic.bphi_fraction
                )
            if np.isfinite(diagnostic.rhs_norm):
                self.maximum_rhs_norm = max(
                    self.maximum_rhs_norm, diagnostic.rhs_norm
                )
            raise
        self.minimum_abs_B_phi = min(self.minimum_abs_B_phi, absolute)
        self.minimum_bphi_fraction = min(self.minimum_bphi_fraction, fraction)
        self.maximum_rhs_norm = max(self.maximum_rhs_norm, norm)
        return rhs


class ToroidalTracer:
    """Trace magnetic field lines using physical toroidal angle ``phi``."""

    def __init__(
        self,
        magnetic_field: MagneticField,
        nfp: int,
        *,
        options: ToroidalTracingOptions | None = None,
    ):
        nfp = operator.index(nfp)
        if nfp < 1:
            raise ValueError("nfp must be positive")
        if not callable(magnetic_field):
            raise TypeError("magnetic_field must be callable")
        self.magnetic_field = magnetic_field
        self.nfp = nfp
        self.options = options or ToroidalTracingOptions()

    @classmethod
    def from_solution(
        cls,
        solution,
        *,
        theta_stride=1,
        zeta_stride=1,
        options=None,
    ):
        """Construct a tracer using a cached REGCOIL field evaluator."""
        magnetic_field = solution.prepare_magnetic_field(
            theta_stride=theta_stride, zeta_stride=zeta_stride
        )
        return cls(magnetic_field, solution.problem.nfp, options=options)

    @property
    def field_period(self):
        return 2 * np.pi / self.nfp

    def rhs(self, phi, rz, *, options=None):
        """Return ``d(R,Z)/dphi`` with conditioning guards applied."""
        return _evaluate_toroidal_rhs(
            phi, rz, self.magnetic_field, options or self.options
        )[0]

    def diagnose(self, phi, rz, *, options=None):
        """Return a complete diagnostic without propagating an invalid point."""
        options = options or self.options
        R, Z = np.asarray(rz, dtype=float)
        try:
            rhs, _, fraction, norm = _evaluate_toroidal_rhs(
                phi, (R, Z), self.magnetic_field, options
            )
        except ToroidalTracingError as error:
            return error.diagnostic
        B_R, B_phi, B_Z = cylindrical_field(
            phi, (R, Z), self.magnetic_field
        )
        magnitude = float(np.sqrt(B_R * B_R + B_phi * B_phi + B_Z * B_Z))
        return ToroidalDiagnostic(
            True,
            "ok",
            float(phi),
            float(R),
            float(Z),
            B_R,
            B_phi,
            B_Z,
            magnitude,
            fraction,
            norm,
            rhs,
        )

    def _resolve_direction(self, diagnostics, direction):
        if direction != "auto":
            value = float(direction)
            if value not in (-1.0, 1.0):
                raise ValueError("direction must be 'auto', -1, or +1")
            return int(value)
        signs = {
            int(np.sign(item.B_phi))
            for item in diagnostics
            if np.isfinite(item.B_phi) and item.B_phi != 0
        }
        if len(signs) > 1:
            raise ValueError(
                "automatic toroidal direction is ambiguous because initial "
                "B_phi values have mixed signs"
            )
        return signs.pop() if signs else 1

    @staticmethod
    def _sampling_schedule(periods, phase_fractions, direction, field_period):
        phase_fractions = np.asarray(phase_fractions, dtype=float)
        if phase_fractions.ndim != 1 or not len(phase_fractions):
            raise ValueError("phase_fractions must be a non-empty 1D sequence")
        if not np.all(np.isfinite(phase_fractions)):
            raise ValueError("phase_fractions must be finite")
        phase_fractions = np.mod(phase_fractions, 1.0)
        if len(np.unique(np.round(phase_fractions, 14))) != len(phase_fractions):
            raise ValueError("phase_fractions must be distinct modulo one period")
        zero = np.flatnonzero(np.isclose(phase_fractions, 0.0))
        if len(zero) != 1:
            raise ValueError("phase_fractions must contain 0 exactly once")

        integration_fractions = np.mod(direction * phase_fractions, 1.0)
        order = np.argsort(integration_fractions)
        integration_fractions = integration_fractions[order]
        ordered_phase_indices = np.arange(len(phase_fractions))[order]

        u = []
        phase_ids = []
        for period in range(periods):
            u.extend(period + integration_fractions)
            phase_ids.extend(ordered_phase_indices)
        u.append(float(periods))
        phase_ids.append(int(zero[0]))
        u = np.asarray(u)
        phase_ids = np.asarray(phase_ids)
        phi = direction * field_period * u
        sample_indices = tuple(
            np.flatnonzero(phase_ids == index)
            for index in range(len(phase_fractions))
        )
        return phase_fractions, phi, sample_indices

    def trace_poincare(
        self,
        initial_rz,
        *,
        periods,
        phase_fractions=(0.0,),
        direction="auto",
        max_workers=1,
        options=None,
    ):
        """Trace independent lines sampled at fixed physical sections.

        Expected conditioning failures are returned as unsuccessful traces.
        Unexpected exceptions propagate so implementation errors are not
        silently mislabeled as bad field lines.
        """
        periods = operator.index(periods)
        if periods < 1:
            raise ValueError("periods must be positive")
        initial_rz = np.asarray(initial_rz, dtype=float)
        if initial_rz.ndim == 1:
            initial_rz = initial_rz[None, :]
        if initial_rz.ndim != 2 or initial_rz.shape[1] != 2:
            raise ValueError("initial_rz must have shape (nlines, 2)")
        if not len(initial_rz):
            raise ValueError("at least one initial condition is required")
        if not np.all(np.isfinite(initial_rz)):
            raise ValueError("initial_rz must contain only finite values")
        options = options or self.options

        initial_diagnostics = tuple(
            self.diagnose(0.0, rz, options=options) for rz in initial_rz
        )
        direction = self._resolve_direction(initial_diagnostics, direction)
        phases, sample_phi, sample_indices = self._sampling_schedule(
            periods, phase_fractions, direction, self.field_period
        )

        def trace_one(rz0):
            monitor = _TraceMonitor(self, options)
            try:
                monitor(0.0, rz0)
                result = solve_ivp(
                    monitor,
                    (float(sample_phi[0]), float(sample_phi[-1])),
                    rz0,
                    method=options.method,
                    rtol=options.rtol,
                    atol=options.atol,
                    t_eval=sample_phi,
                    dense_output=False,
                )
            except ToroidalTracingError as error:
                return PoincareTrace(
                    np.array(rz0, copy=True),
                    False,
                    np.empty(0),
                    np.empty((0, 2)),
                    0,
                    str(error),
                    error.diagnostic,
                    float(monitor.minimum_abs_B_phi),
                    float(monitor.minimum_bphi_fraction),
                    float(monitor.maximum_rhs_norm),
                    phases,
                    sample_indices,
                    direction,
                    self.nfp,
                    periods,
                )

            return PoincareTrace(
                np.array(rz0, copy=True),
                bool(result.success),
                np.asarray(result.t),
                np.asarray(result.y).T,
                int(result.nfev),
                "ok" if result.success else str(result.message),
                monitor.last_diagnostic,
                float(monitor.minimum_abs_B_phi),
                float(monitor.minimum_bphi_fraction),
                float(monitor.maximum_rhs_norm),
                phases,
                sample_indices,
                direction,
                self.nfp,
                periods,
            )

        max_workers = operator.index(max_workers)
        if max_workers < 1:
            raise ValueError("max_workers must be positive")
        start = perf_counter()
        if max_workers == 1 or len(initial_rz) == 1:
            traces = tuple(trace_one(rz) for rz in initial_rz)
        else:
            with ThreadPoolExecutor(max_workers=max_workers) as executor:
                traces = tuple(executor.map(trace_one, initial_rz))
        wall_time = perf_counter() - start
        return PoincareBatch(
            traces,
            phases,
            direction,
            self.nfp,
            periods,
            wall_time,
        )

    def first_return_map(self, rz, *, direction=1, options=None):
        """Map one ``phi=0`` point through exactly one field period."""
        if direction not in (-1, 1):
            raise ValueError("direction must be -1 or +1")
        options = options or self.options
        result = solve_ivp(
            lambda phi, state: self.rhs(phi, state, options=options),
            (0.0, direction * self.field_period),
            np.asarray(rz, dtype=float),
            method=options.method,
            rtol=options.rtol,
            atol=options.atol,
            dense_output=False,
        )
        if not result.success:
            raise RuntimeError(result.message)
        return result.y[:, -1]

    def find_magnetic_axis(
        self,
        plasma,
        *,
        direction=1,
        initial_guess=None,
        fixed_point_tolerance=1.0e-8,
        jacobian_relative_step=1.0e-5,
        options=None,
        max_nfev=30,
    ):
        """Find and validate a local elliptic fixed point at physical phi=0."""
        from matplotlib.path import Path

        options = options or ToroidalTracingOptions(rtol=1e-10, atol=1e-11)
        R_section, Z_section = plasma.cross_section(phi=0.0)
        R_boundary = np.asarray(R_section).squeeze()
        Z_boundary = np.asarray(Z_section).squeeze()
        section_path = Path(np.column_stack([R_boundary, Z_boundary]))
        lower = np.array([R_boundary.min(), Z_boundary.min()])
        upper = np.array([R_boundary.max(), Z_boundary.max()])
        scale = upper - lower
        if np.any(scale <= 0):
            raise ValueError("plasma cross-section has a degenerate bounding box")
        if initial_guess is None:
            initial_guess = 0.5 * (lower + upper)
        initial_guess = np.asarray(initial_guess, dtype=float)
        if initial_guess.shape != (2,):
            raise ValueError("initial_guess must have shape (2,)")
        if not section_path.contains_point(initial_guess):
            raise ValueError("initial_guess is outside the plasma cross-section")

        evaluations = 0

        def residual(rz):
            nonlocal evaluations
            evaluations += 1
            return self.first_return_map(
                rz, direction=direction, options=options
            ) - rz

        def residual_jacobian(rz):
            jacobian = np.empty((2, 2))
            for column, step in enumerate(jacobian_relative_step * scale):
                offset = np.zeros(2)
                offset[column] = step
                jacobian[:, column] = (
                    residual(rz + offset) - residual(rz - offset)
                ) / (2 * step)
            return jacobian

        start = perf_counter()
        fit = least_squares(
            residual,
            initial_guess,
            jac=residual_jacobian,
            bounds=(lower, upper),
            x_scale=scale,
            xtol=1.0e-11,
            ftol=1.0e-12,
            gtol=1.0e-12,
            max_nfev=max_nfev,
        )
        axis_rz = fit.x
        return_residual = residual(axis_rz)
        return_error = float(np.linalg.norm(return_residual))
        if not np.isfinite(return_error) or return_error > fixed_point_tolerance:
            raise RuntimeError(
                "fixed-point solve did not meet the requested tolerance: "
                f"return_error={return_error:.3e} m; {fit.message}"
            )
        if not section_path.contains_point(axis_rz):
            raise RuntimeError("fixed-point candidate is outside the plasma")

        return_map_jacobian = residual_jacobian(axis_rz) + np.eye(2)
        eigenvalues = np.linalg.eigvals(return_map_jacobian)
        determinant = float(np.linalg.det(return_map_jacobian))
        discriminant = float(
            np.trace(return_map_jacobian) ** 2 - 4 * determinant
        )
        if determinant <= 0 or discriminant >= 0:
            raise RuntimeError(
                "a fixed point was found, but its return-map linearization "
                "is not elliptic"
            )
        return MagneticAxisResult(
            np.asarray(axis_rz),
            np.asarray(initial_guess),
            return_residual,
            return_error,
            return_map_jacobian,
            eigenvalues,
            determinant,
            discriminant,
            evaluations,
            str(fit.message),
            perf_counter() - start,
        )


def estimate_iota(
    batch: PoincareBatch,
    axis_rz,
    *,
    period_counts: Sequence[int] = (24, 50, 100, 200),
    minimum_radius=1.0e-6,
):
    """Estimate geometric rotational transform from repeated phi=0 returns."""
    axis_rz = np.asarray(axis_rz, dtype=float)
    if axis_rz.shape != (2,):
        raise ValueError("axis_rz must have shape (2,)")
    counts = np.asarray(
        [operator.index(count) for count in period_counts], dtype=int
    )
    counts = counts[(counts > 0) & (counts <= batch.periods)]
    if len(counts) == 0:
        raise ValueError("period_counts contains no values within the batch")

    rows = []
    for trace in batch.successful:
        radius = float(np.linalg.norm(trace.initial_rz - axis_rz))
        if radius < minimum_radius:
            continue
        return_phi, return_rz = trace.same_section_returns()
        angle = np.unwrap(np.arctan2(
            return_rz[:, 1] - axis_rz[1],
            return_rz[:, 0] - axis_rz[0],
        ))
        estimates = np.array([
            np.polyfit(return_phi[: count + 1], angle[: count + 1], 1)[0]
            for count in counts
        ])
        rows.append((radius, trace.initial_rz, estimates))

    rows.sort(key=lambda row: row[0])
    if rows:
        radius = np.array([row[0] for row in rows])
        initial = np.stack([row[1] for row in rows])
        iota = np.stack([row[2] for row in rows])
    else:
        radius = np.empty(0)
        initial = np.empty((0, 2))
        iota = np.empty((0, len(counts)))
    return IotaProfile(radius, initial, counts, iota, axis_rz)


__all__ = [
    "IotaProfile",
    "MagneticAxisResult",
    "PoincareBatch",
    "PoincareTrace",
    "ToroidalDiagnostic",
    "ToroidalTracer",
    "ToroidalTracingError",
    "ToroidalTracingOptions",
    "cylindrical_field",
    "estimate_iota",
]
