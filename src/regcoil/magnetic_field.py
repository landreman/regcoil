"""Prepared direct-quadrature magnetic-field evaluation.

The expensive, position-independent winding-surface data are assembled once.
Singleton calls use the fused C kernel used by adaptive field-line integrators;
batched calls use bounded NumPy temporaries.
"""

from __future__ import annotations

from dataclasses import dataclass, field
import operator
from typing import TYPE_CHECKING

import numpy as np

from . import _core

if TYPE_CHECKING:
    from .regcoil import Solution


DEFAULT_CHUNK_SIZE = 64


def _flatten_grid(array):
    """Flatten ``(..., ntheta, nzeta)`` with theta varying fastest."""
    return np.moveaxis(array, -2, -1).reshape(*array.shape[:-2], -1)


@dataclass(frozen=True, eq=False)
class MagneticFieldEvaluator:
    """Prepared direct-quadrature evaluator for one REGCOIL solution.

    Instances are normally obtained from
    :meth:`regcoil.Solution.prepare_magnetic_field` rather than constructed
    directly.
    """

    source_points: np.ndarray = field(repr=False)
    source_radius_squared: np.ndarray = field(repr=False)
    weighted_surface_current: np.ndarray = field(repr=False)
    weighted_current_moment: np.ndarray = field(repr=False)
    _source_points_c: np.ndarray = field(init=False, repr=False)
    _weighted_surface_current_c: np.ndarray = field(init=False, repr=False)
    theta_stride: int = 1
    zeta_stride: int = 1
    quadrature_shape: tuple[int, int] = field(default=(0, 0))

    def __post_init__(self):
        # The compiled singleton kernel requires row-major (nsource, 3)
        # arrays. Normally these are views; direct construction may copy once.
        object.__setattr__(
            self,
            "_source_points_c",
            np.ascontiguousarray(self.source_points, dtype=np.float64),
        )
        object.__setattr__(
            self,
            "_weighted_surface_current_c",
            np.ascontiguousarray(
                self.weighted_surface_current, dtype=np.float64
            ),
        )

    @property
    def nsource(self):
        """Number of retained full-torus quadrature points."""
        return self.source_points.shape[0]

    @staticmethod
    def _validate_stride(name, stride, grid_size):
        if isinstance(stride, bool):
            raise TypeError(f"{name} must be an integer, got {stride!r}")
        try:
            stride = operator.index(stride)
        except TypeError as exc:
            raise TypeError(
                f"{name} must be an integer, got {stride!r}"
            ) from exc
        if stride < 1:
            raise ValueError(f"{name} must be positive, got {stride}")
        if grid_size % stride:
            raise ValueError(
                f"{name}={stride} must divide its periodic grid size "
                f"{grid_size}"
            )
        return stride

    @classmethod
    def _validated_strides(cls, solution, theta_stride, zeta_stride):
        coil = solution.problem.coil
        return (
            cls._validate_stride(
                "theta_stride", theta_stride, coil.ntheta
            ),
            cls._validate_stride(
                "zeta_stride", zeta_stride, coil.nfp * coil.nzeta
            ),
        )

    @classmethod
    def from_solution(
        cls,
        solution: Solution,
        *,
        theta_stride=1,
        zeta_stride=1,
    ):
        """Prepare a field evaluator from one current-potential solution."""
        coil = solution.problem.coil
        current_one_period = solution.current_density()
        theta_stride, zeta_stride = cls._validated_strides(
            solution, theta_stride, zeta_stride
        )

        source_blocks = []
        weighted_current_blocks = []
        weights_one_period = (
            _flatten_grid(coil.norm_normal) * coil.dtheta * coil.dzeta
        )

        # coil.r spans the full torus, while current_density() covers one
        # field period. Generate the remaining currents by rigid rotation.
        for period in range(coil.nfp):
            start = period * coil.nzeta
            stop = start + coil.nzeta
            source_period = _flatten_grid(coil.r[:, :, start:stop]).T

            angle = 2 * np.pi * period / coil.nfp
            cosine = np.cos(angle)
            sine = np.sin(angle)
            current_period = np.empty_like(current_one_period)
            current_period[0] = (
                cosine * current_one_period[0]
                - sine * current_one_period[1]
            )
            current_period[1] = (
                sine * current_one_period[0]
                + cosine * current_one_period[1]
            )
            current_period[2] = current_one_period[2]
            current_period = _flatten_grid(current_period).T

            source_blocks.append(source_period)
            weighted_current_blocks.append(
                current_period * weights_one_period[:, None]
            )

        full_grid_shape = (coil.nfp * coil.nzeta, coil.ntheta, 3)
        source_points = np.concatenate(source_blocks, axis=0).reshape(
            full_grid_shape
        )[::zeta_stride, ::theta_stride].reshape(-1, 3)
        weighted_surface_current = np.concatenate(
            weighted_current_blocks, axis=0
        ).reshape(full_grid_shape)[
            ::zeta_stride, ::theta_stride
        ].reshape(-1, 3)
        weighted_surface_current = (
            weighted_surface_current * theta_stride * zeta_stride
        )

        return cls(
            source_points=source_points,
            source_radius_squared=np.einsum(
                "ij,ij->i", source_points, source_points
            ),
            weighted_surface_current=weighted_surface_current,
            weighted_current_moment=np.cross(
                weighted_surface_current, source_points
            ),
            theta_stride=theta_stride,
            zeta_stride=zeta_stride,
            quadrature_shape=(
                coil.nfp * coil.nzeta // zeta_stride,
                coil.ntheta // theta_stride,
            ),
        )

    def __call__(self, points, *, chunk_size=DEFAULT_CHUNK_SIZE):
        """Evaluate the surface-current magnetic field at ``points``.

        ``chunk_size`` bounds the main temporary to approximately
        ``8 * chunk_size * nsource`` bytes. Direct quadrature is singular on,
        and is not intended for evaluation very close to, the winding surface.
        """
        if isinstance(chunk_size, bool):
            raise TypeError(
                f"chunk_size must be an integer, got {chunk_size!r}"
            )
        try:
            chunk_size = operator.index(chunk_size)
        except TypeError as exc:
            raise TypeError(
                f"chunk_size must be an integer, got {chunk_size!r}"
            ) from exc
        if chunk_size < 1:
            raise ValueError("chunk_size must be positive")

        points = np.asarray(points, dtype=np.float64)
        if points.ndim == 0 or points.shape[-1] != 3:
            raise ValueError("points must have shape (..., 3)")
        if not np.all(np.isfinite(points)):
            raise ValueError("points must contain only finite values")

        original_shape = points.shape
        points = points.reshape(-1, 3)
        if len(points) == 1:
            return self._evaluate_single(points[0]).reshape(original_shape)

        magnetic_field = np.empty_like(points)
        for start in range(0, len(points), chunk_size):
            stop = min(start + chunk_size, len(points))
            point_block = points[start:stop]
            if len(point_block) == 1:
                magnetic_field[start] = self._evaluate_single(point_block[0])
                continue

            distance_squared = -2.0 * (point_block @ self.source_points.T)
            distance_squared += np.einsum(
                "ij,ij->i", point_block, point_block
            )[:, None]
            distance_squared += self.source_radius_squared[None, :]
            if np.any(distance_squared <= 0.0):
                raise ValueError(
                    "Biot-Savart evaluation point lies on a coil-surface "
                    "quadrature point"
                )

            np.power(distance_squared, -1.5, out=distance_squared)
            summed_current = distance_squared @ self.weighted_surface_current
            summed_moment = distance_squared @ self.weighted_current_moment
            magnetic_field[start:stop] = 1.0e-7 * (
                np.cross(summed_current, point_block) - summed_moment
            )

        return magnetic_field.reshape(original_shape)

    def _evaluate_single(self, point):
        """Fast path for singleton calls made by ODE integrators."""
        return _core.magnetic_field_single(  # type: ignore[attr-defined]
            self._source_points_c,
            self._weighted_surface_current_c,
            np.ascontiguousarray(point),
        )
