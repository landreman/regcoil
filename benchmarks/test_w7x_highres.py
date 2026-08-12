"""CodSpeed benchmarks for the quickstart W7-X problem at high resolution.

This is the README/quickstart case with ntheta = nzeta = 128 on both surfaces,
mpol = ntor = 32 for the winding surface, and mpol_potential = ntor_potential
= 32 for the current potential.  Three steps are timed separately:

* ``CoilSurface.from_uniform_offset()`` -- building the winding surface,
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
NTHETA = 128
NZETA = 128
MPOL = 32
NTOR = 32
MPOL_POTENTIAL = 32
NTOR_POTENTIAL = 32
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


def _make_problem(plasma, coil):
    return regcoil.Regcoil(
        plasma, coil,
        mpol_potential=MPOL_POTENTIAL, ntor_potential=NTOR_POTENTIAL,
    )


@pytest.fixture(scope="session")
def coil(plasma):
    """Winding surface, built once and shared by the later benchmarks."""
    return _make_coil(plasma)


@pytest.fixture(scope="session")
def problem(plasma, coil):
    """Regcoil object with its operators already assembled."""
    return _make_problem(plasma, coil)


def test_from_uniform_offset(benchmark, plasma):
    coil = benchmark(_make_coil, plasma)
    assert coil.ntheta == NTHETA


def test_regcoil_init(benchmark, plasma, coil):
    problem = benchmark(_make_problem, plasma, coil)
    assert problem.nbf > 0


def test_solve(benchmark, problem):
    solution = benchmark(problem.solve, lam=LAM)
    assert solution.f_B > 0
