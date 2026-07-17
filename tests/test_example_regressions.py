"""Run example regression checks via pytest (mirrors ``make test``).

Copies of each ``examples/*/tests.py`` live under ``tests/regression/``.
Inputs remain in ``examples/`` (namelists use ``../../equilibria/...`` paths).
"""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
EXAMPLES_DIR = REPO_ROOT / "examples"
REGRESSION_DIR = Path(__file__).resolve().parent / "regression"
REGCOIL_EXE = REPO_ROOT / "regcoil"


def _example_names() -> list[str]:
    names: list[str] = []
    for tests_py in sorted(REGRESSION_DIR.glob("*/tests.py")):
        name = tests_py.parent.name
        nml = EXAMPLES_DIR / name / f"regcoil_in.{name}"
        if nml.is_file():
            names.append(name)
    return names


EXAMPLE_NAMES = _example_names()


def _ntheta_plasma(example_name: str) -> int | None:
    """Return the active ``ntheta_plasma`` from the example namelist, if set."""
    nml = EXAMPLES_DIR / example_name / f"regcoil_in.{example_name}"
    for raw in nml.read_text().splitlines():
        line = raw.split("!")[0].strip()
        if not line.lower().startswith("ntheta_plasma"):
            continue
        # e.g. ntheta_plasma = 128  or  ntheta_plasma=128
        rhs = line.split("=", 1)[1].strip().rstrip(",").strip()
        return int(float(rhs))
    return None


@pytest.fixture(scope="session")
def regcoil_executable() -> Path:
    if not REGCOIL_EXE.is_file() or not os.access(REGCOIL_EXE, os.X_OK):
        pytest.skip(
            "regcoil executable not found; build with `make` before running "
            "example regression tests"
        )
    return REGCOIL_EXE


@pytest.mark.parametrize("example_name", EXAMPLE_NAMES)
def test_example_regression(example_name: str, regcoil_executable: Path) -> None:
    example_dir = EXAMPLES_DIR / example_name
    assert example_dir.is_dir(), f"missing example directory {example_dir}"

    ntheta = _ntheta_plasma(example_name)
    if ntheta is not None and ntheta >= 128:
        pytest.skip(f"high-resolution example (ntheta_plasma={ntheta})")

    output_nc = example_dir / f"regcoil_out.{example_name}.nc"
    if output_nc.exists():
        output_nc.unlink()

    input_file = f"regcoil_in.{example_name}"
    submit = os.environ.get("REGCOIL_COMMAND_TO_SUBMIT_JOB", "").split()
    cmd = submit + [str(regcoil_executable), input_file]
    run = subprocess.run(
        cmd,
        cwd=example_dir,
        capture_output=True,
        text=True,
    )
    if run.returncode != 0:
        pytest.fail(
            f"regcoil failed for {example_name} (exit {run.returncode})\n"
            f"command: {cmd}\n"
            f"stdout:\n{run.stdout}\n"
            f"stderr:\n{run.stderr}"
        )

    tests_py = REGRESSION_DIR / example_name / "tests.py"
    check = subprocess.run(
        [sys.executable, str(tests_py)],
        cwd=example_dir,
        capture_output=True,
        text=True,
    )
    if check.returncode != 0:
        pytest.fail(
            f"regression checks failed for {example_name} "
            f"(exit {check.returncode})\n"
            f"stdout:\n{check.stdout}\n"
            f"stderr:\n{check.stderr}"
        )
