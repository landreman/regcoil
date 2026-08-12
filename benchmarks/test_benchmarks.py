"""CodSpeed benchmarks for the quickstart W7-X problem at high resolution.

This is the README/quickstart case at the resolution set by the constants
below, with four steps timed separately:

* ``CoilSurface.from_uniform_offset()`` -- building the winding surface,
* the same call with ``standard_toroidal_angle=True``, the legacy
  root-solve construction,
* ``Regcoil()`` -- assembling the inductance / regularization operators,
* ``Regcoil.solve()`` -- the linear solve and diagnostics for one lambda.

These files are deliberately outside ``tests/`` (``testpaths`` in
pyproject.toml) so a plain ``pytest`` run never picks up a multi-minute
benchmark.  Run them by hand with::

    pip install pytest-codspeed
    pytest benchmarks/ --codspeed

In CI they are run by .github/workflows/codspeed.yml, which uploads the
walltimes to CodSpeed so a pull request is compared against master.
"""

import pytest

import regcoil

# Resolution of the benchmark problem.  Bumping any of these invalidates the
# historical CodSpeed series for these benchmarks, so change them deliberately.
NTHETA = 64
NZETA = 64
MPOL = 12
NTOR = 12
MPOL_POTENTIAL = 12
NTOR_POTENTIAL = 12
SEPARATION = 0.3
LAM = 1e-14


@pytest.fixture(scope="session")
def plasma():
    """The W7-X plasma boundary with virtual-casing B_normal."""
    ds = regcoil.examples("W7-X")
    surface = regcoil.PlasmaSurface.from_wout(ds.wout, ntheta=NTHETA, nzeta=NZETA)
    surface.set_bnormal_from_virtual_casing(ds.vcasing)
    return surface


def _make_coil(plasma):
    return regcoil.CoilSurface.from_uniform_offset(
        plasma, separation=SEPARATION, ntheta=NTHETA, nzeta=NZETA,
        mpol=MPOL, ntor=NTOR,
    )


def _make_coil_standard_angle(plasma):
    """Legacy `geometry_option_coil=2` construction: every grid point is a
    root solve for the offset-surface point whose Cartesian toroidal angle
    equals the coil `zeta`, rather than a direct move along the normal."""
    return regcoil.CoilSurface.from_uniform_offset(
        plasma, separation=SEPARATION, ntheta=NTHETA, nzeta=NZETA,
        mpol=MPOL, ntor=NTOR, standard_toroidal_angle=True,
    )


def _make_problem(plasma, coil):
    return regcoil.Regcoil(
        plasma, coil,
        mpol_potential=MPOL_POTENTIAL, ntor_potential=NTOR_POTENTIAL,
    )


@pytest.fixture(scope="session")
def coil(plasma):
    """Winding surface used as input to the later benchmarks."""
    return _make_coil(plasma)


@pytest.fixture(scope="session")
def problem(plasma, coil):
    """Regcoil object with its operators already assembled."""
    return _make_problem(plasma, coil)


def test_from_uniform_offset(benchmark, plasma):
    coil = benchmark(_make_coil, plasma)
    assert coil.ntheta == NTHETA


def test_offset_standard_zeta(benchmark, plasma):
    coil = benchmark(_make_coil_standard_angle, plasma)
    assert coil.ntheta == NTHETA


def test_regcoil_init(benchmark, plasma, coil):
    problem = benchmark(_make_problem, plasma, coil)
    assert problem.nbf > 0


def test_solve(benchmark, problem):
    solution = benchmark(problem.solve, lam=LAM)
    assert solution.f_B > 0
