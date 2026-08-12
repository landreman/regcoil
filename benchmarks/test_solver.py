"""Benchmarks for the solver itself.

`test_operator_assembly` is the one that dominates a production run: the
Fortran `build_g_and_h` kernel, the `matrix_B`/`matrix_K` BLAS level-3 calls,
and the generalized eigensolve.  Everything after it is O(nbf**2) per lambda
thanks to the cached eigendecomposition, and the benchmarks below keep that
property honest.
"""

from __future__ import annotations

import numpy as np

from ._problem import LAM, LAMBDAS, make_problem


def test_operator_assembly(benchmark, plasma, coil):
    """`Regcoil(...)`: Fortran kernel + matrix assembly + generalized eigensolve."""
    problem = benchmark(make_problem, plasma, coil)
    assert problem.nbf > 0


def test_solve_single_lambda(benchmark, problem):
    """One regularized solve on top of the cached eigendecomposition, including
    the `Solution` diagnostics (f_B, f_K, max_K, ...)."""
    solution = benchmark(problem.solve, lam=LAM)
    assert solution.f_B > 0


def test_solve_infinite_lambda(benchmark, problem):
    """The heavily-regularized limit, which takes the separate `lam=inf`
    branch of the coefficient formula."""
    solution = benchmark(problem.solve, lam=np.inf)
    assert solution.f_K > 0


def test_lambda_scan(benchmark, problem):
    """A whole regularization scan: vectorized coefficients, then the
    per-lambda diagnostics."""
    scan = benchmark(problem.scan, LAMBDAS)
    assert len(scan) == LAMBDAS.size


def test_solve_for_target(benchmark, problem):
    """Bisection in log(lambda) for a prescribed `max_K` -- a handful of solves
    driven by `scipy.optimize.brentq`."""
    lo = problem.solve(0.0).max_K
    hi = problem.solve(np.inf).max_K
    target = 0.5 * (lo + hi)
    solution = benchmark(problem.solve_for_target, "max_K", target)
    assert np.isclose(solution.max_K, target, rtol=1e-6)


def test_current_potential(benchmark, solution):
    """Expand the current potential from mode amplitudes onto the coil grid."""
    Phi = benchmark(solution.current_potential)
    assert Phi.shape == (solution.problem.coil.ntheta, solution.problem.coil.nzeta)


def test_current_density(benchmark, solution):
    """Expand the Cartesian surface current density onto the coil grid (a
    contraction against the full `f_all` operator)."""
    K = benchmark(solution.current_density)
    assert K.shape == (3, solution.problem.coil.ntheta, solution.problem.coil.nzeta)
