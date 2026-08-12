"""Benchmarks for the geometry front end: reading boundaries, evaluating a
Fourier surface, and building a winding surface.

`CoilSurface.from_uniform_offset` is the expensive one -- it samples the offset
point cloud, builds the theta map, then Fourier-transforms the result -- and it
is pure numpy, so it is also the part most likely to drift when the geometry
code is refactored.
"""

from __future__ import annotations

import numpy as np
import pytest

import regcoil

from ._problem import (
    MPOL,
    NTHETA,
    NTHETA_COARSE,
    NTOR,
    NZETA,
    NZETA_COARSE,
    SEPARATION,
    make_coil,
    make_plasma,
)


def test_plasma_from_wout(benchmark, w7x):
    """Read the VMEC wout file and build the plasma boundary (which also
    evaluates |B| on the boundary grid)."""
    plasma = benchmark(
        regcoil.PlasmaSurface.from_wout, w7x.wout, ntheta=NTHETA, nzeta=NZETA
    )
    assert plasma.nfp == 5


def test_plasma_from_wout_with_bnormal(benchmark, w7x):
    """The full quickstart plasma-side setup: boundary plus the virtual-casing
    B_normal read from file."""
    plasma = benchmark(make_plasma, w7x)
    assert plasma.Bnormal_from_plasma_current.shape == (NTHETA, NZETA)


def test_plasma_bnormal_from_bnorm_file(benchmark, w7x):
    """The stellopt BNORM route to `Bnormal_from_plasma_current`: an ASCII mode
    list summed onto the plasma grid."""
    plasma = regcoil.PlasmaSurface.from_wout(w7x.wout, ntheta=NTHETA, nzeta=NZETA)
    benchmark(plasma.set_bnormal_from_bnorm_file, w7x.bnorm)
    assert np.any(plasma.Bnormal_from_plasma_current != 0)


def test_surface_evaluate_at(benchmark, plasma):
    """Evaluation at arbitrary paired (theta, zeta) points -- the inner loop of
    coil cutting and of every reparameterization."""
    theta_2d, zeta_2d = np.meshgrid(plasma.theta, plasma.zeta, indexing="ij")
    theta_pts = np.ravel(theta_2d)
    zeta_pts = np.ravel(zeta_2d)
    data = benchmark(plasma.evaluate_at, theta_pts, zeta_pts)
    assert data["r"].shape == (3, theta_pts.size)


def test_coil_from_nescin(benchmark, w7x, plasma):
    """Read a pre-computed winding surface in nescin format."""
    coil = benchmark(
        regcoil.CoilSurface.from_nescin,
        w7x.nescin, nfp=plasma.nfp, ntheta=NTHETA, nzeta=NZETA,
    )
    assert coil.nfp == plasma.nfp


def test_coil_from_uniform_offset(benchmark, plasma):
    """Uniform offset along the plasma normal, with the default
    uniform-arclength theta reparameterization."""
    coil = benchmark(make_coil, plasma)
    assert coil.ntheta == NTHETA


def test_coil_from_uniform_offset_no_reparameterization(benchmark, plasma):
    """The same construction with the theta map switched off, so the difference
    from the benchmark above is the cost of the reparameterization."""
    coil = benchmark(
        regcoil.CoilSurface.from_uniform_offset,
        plasma, separation=SEPARATION, ntheta=NTHETA, nzeta=NZETA,
        mpol=MPOL, ntor=NTOR, theta_reparameterization=None,
    )
    assert coil.theta_map is None


def test_coil_from_uniform_offset_standard_toroidal_angle(benchmark, plasma):
    """Legacy `geometry_option_coil=4`: a scalar root solve per grid point to
    place each offset point at a prescribed toroidal angle.  Run on a coarser
    grid than the other offset benchmarks because it costs several times more
    per point.
    """
    coil = benchmark(
        regcoil.CoilSurface.from_uniform_offset,
        plasma, separation=SEPARATION,
        ntheta=NTHETA_COARSE, nzeta=NZETA_COARSE, mpol=MPOL, ntor=NTOR,
        standard_toroidal_angle=True,
    )
    assert coil.standard_toroidal_angle


@pytest.mark.parametrize(
    "scheme",
    [
        pytest.param(regcoil.UniformArclength(), id="uniform_arclength"),
        pytest.param(regcoil.CurvatureWeighted(), id="curvature_weighted"),
    ],
)
def test_reparameterize_theta(benchmark, plasma, scheme):
    """Poloidal-angle reparameterization of an existing surface: quadrature for
    the theta map, resampling, then a refit."""
    surface = benchmark(
        plasma.reparameterize_theta,
        scheme, mpol=MPOL, ntor=NTOR, ntheta=NTHETA, nzeta=NZETA,
    )
    assert surface.theta_map is not None
