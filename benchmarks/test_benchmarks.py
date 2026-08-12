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
MPOL = 18
NTOR = 18
MPOL_POTENTIAL = 18
NTOR_POTENTIAL = 18
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


# Objects handed down from one benchmark to the next.  `benchmark(f, ...)`
# returns f's real return value, so the coil surface built while benchmarking
# `from_uniform_offset` is the same object the `Regcoil` benchmark needs, and
# likewise for the Regcoil object and `solve`.  Reusing them keeps the wall
# clock of the whole job down: a fresh `Regcoil()` here would cost as much as
# a whole extra measured round.  The fixtures below fall back to building the
# object if its producing benchmark did not run (e.g. `pytest -k solve`).
_built: dict = {}


@pytest.fixture(scope="session")
def coil(plasma):
    """Winding surface, from test_from_uniform_offset if it has run."""
    if "coil" not in _built:
        _built["coil"] = _make_coil(plasma)
    return _built["coil"]


@pytest.fixture(scope="session")
def problem(plasma, coil):
    """Regcoil object with its operators assembled, from test_regcoil_init."""
    if "problem" not in _built:
        _built["problem"] = _make_problem(plasma, coil)
    return _built["problem"]


# NOTE ON CALL COUNTS: pytest-codspeed's walltime instrument calls the measured
# function at least three times -- once untimed to capture the return value,
# at least once for warmup, then at least one timed round.  Only that last
# group appears in the "Run time" column of the results table, so a benchmark
# reported as 141 s costs upwards of 7 minutes of job time.  The floor of three
# is structural in the plugin and cannot be configured away.


def test_from_uniform_offset(benchmark, plasma):
    _built["coil"] = benchmark(_make_coil, plasma)
    assert _built["coil"].ntheta == NTHETA


def test_from_uniform_offset_standard_toroidal_angle(benchmark, plasma):
    # Not fed into the chain below: the rest of the benchmarks deliberately
    # stay on the default (normal-offset) winding surface.
    coil = benchmark(_make_coil_standard_angle, plasma)
    assert coil.ntheta == NTHETA


def test_regcoil_init(benchmark, plasma, coil):
    _built["problem"] = benchmark(_make_problem, plasma, coil)
    assert _built["problem"].nbf > 0


def test_solve(benchmark, problem):
    solution = benchmark(problem.solve, lam=LAM)
    assert solution.f_B > 0
