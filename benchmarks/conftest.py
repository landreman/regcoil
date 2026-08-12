"""Session-scoped fixtures shared by the benchmark modules.

The expensive setup (plasma surface, winding surface, assembled operators, one
solve) is built once per session so that each benchmark measures only its own
step.  See `benchmarks/_problem.py` for the problem definition and its size.
"""

from __future__ import annotations

import pytest

import regcoil

from ._problem import LAM, make_coil, make_plasma, make_problem


@pytest.fixture(scope="session")
def w7x():
    """Paths to the bundled W7-X example dataset."""
    return regcoil.examples("W7-X")


@pytest.fixture(scope="session")
def plasma(w7x):
    """Plasma surface, built once and shared by every benchmark.

    Benchmarks must not mutate it: its `Bnormal_from_plasma_current` feeds the
    `problem` fixture below.
    """
    return make_plasma(w7x)


@pytest.fixture(scope="session")
def coil(plasma):
    """Winding surface, built once and shared by the solver benchmarks."""
    return make_coil(plasma)


@pytest.fixture(scope="session")
def problem(plasma, coil):
    """`Regcoil` with its operators already assembled."""
    return make_problem(plasma, coil)


@pytest.fixture(scope="session")
def solution(problem):
    """One solve at `LAM`, reused by the post-processing benchmarks."""
    return problem.solve(lam=LAM)
