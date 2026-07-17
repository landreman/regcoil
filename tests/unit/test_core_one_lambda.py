"""Phase 4: one-λ solve via the Fortran extension matches a known example."""

from pathlib import Path

import pytest

import regcoil._core as core

EXAMPLE_DIR = (
    Path(__file__).resolve().parents[2]
    / "examples"
    / "axisymmetrySanityTest_chi2K_regularization"
)
NML = EXAMPLE_DIR / "regcoil_in.axisymmetrySanityTest_chi2K_regularization"

# Same absolute tolerance as examples/axisymmetrySanityTest_chi2K_regularization/tests.py
ABS_TOL = 1e-10


@pytest.fixture(scope="module")
def axisymmetry_setup():
    if not NML.is_file():
        pytest.skip(f"missing example namelist: {NML}")
    core.set_verbose(False)
    core.setup(str(NML))
    return core


def test_import_core():
    assert hasattr(core, "setup")
    assert hasattr(core, "solve_lambda")
    assert hasattr(core, "build_matrices")
    assert hasattr(core, "prepare_solve")


def test_one_lambda_axisymmetry(axisymmetry_setup):
    """Circular plasma/coil surfaces: chi2_B and max_Bnormal vanish at one λ."""
    chi2_B, chi2_K, max_Bnormal, max_K = axisymmetry_setup.solve_lambda(1e-15)
    assert abs(chi2_B) < ABS_TOL
    assert abs(max_Bnormal) < ABS_TOL
    assert chi2_K >= 0.0
    assert max_K >= 0.0


def test_solve_ilambda_axisymmetry(axisymmetry_setup):
    """Namelist lambda index 0 also yields vanishing B-normal metrics."""
    chi2_B, _chi2_K, max_Bnormal, _max_K = axisymmetry_setup.solve_ilambda(0)
    assert abs(chi2_B) < ABS_TOL
    assert abs(max_Bnormal) < ABS_TOL
