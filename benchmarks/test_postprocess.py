"""Benchmarks for what happens after a solve: cutting discrete coils out of the
current potential, and the NetCDF save/load round trip.
"""

from __future__ import annotations

import pytest

import regcoil

from ._problem import COILS_PER_HALF_PERIOD, THICKNESS


def test_cut_coils(benchmark, solution):
    """Trace the current-potential contours and evaluate the coil curves in
    Cartesian space (contourpy plus `Surface.evaluate_at`)."""
    coils = benchmark(solution.cut, COILS_PER_HALF_PERIOD)
    assert len(coils) == 2 * COILS_PER_HALF_PERIOD * solution.problem.nfp


def test_cut_coils_with_thickness(benchmark, solution):
    """The same cut plus the finite-thickness ribbon geometry for every coil."""
    coils = benchmark(solution.cut, COILS_PER_HALF_PERIOD, thickness=THICKNESS)
    assert len(coils.ribbons()) == len(coils)


@pytest.fixture(scope="session")
def saved_solution(solution, tmp_path_factory):
    """A solution written to disk, for the load benchmark."""
    path = tmp_path_factory.mktemp("regcoil_bench") / "solution.nc"
    solution.save(path)
    return path


def test_save_solution(benchmark, solution, tmp_path_factory):
    """Serialize the plasma, coil, problem, and solution to NetCDF-4."""
    path = tmp_path_factory.mktemp("regcoil_bench_save") / "solution.nc"
    benchmark(solution.save, path)
    assert path.exists()


def test_load_solution(benchmark, saved_solution):
    """Read the file back: surfaces plus the stored solution grids, without
    rebuilding the operators."""
    loaded = benchmark(regcoil.load, saved_solution)
    assert loaded.solutions is not None
