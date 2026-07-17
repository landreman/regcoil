from pathlib import Path

import pytest

from regcoil import RegcoilProblem

EXAMPLE_CHI2K = (
    Path(__file__).resolve().parents[2]
    / "examples"
    / "axisymmetrySanityTest_chi2K_regularization"
)
NML_CHI2K = EXAMPLE_CHI2K / "regcoil_in.axisymmetrySanityTest_chi2K_regularization"

EXAMPLE_LB = (
    Path(__file__).resolve().parents[2]
    / "examples"
    / "axisymmetrySanityTest_Laplace_Beltrami_regularization"
)
NML_LB = EXAMPLE_LB / "regcoil_in.axisymmetrySanityTest_Laplace_Beltrami_regularization"

ABS_TOL = 1e-10


@pytest.fixture
def chi2k_problem():
    if not NML_CHI2K.is_file():
        pytest.skip(f"missing example namelist: {NML_CHI2K}")
    prob = RegcoilProblem()
    prob.set_verbose(False)
    prob.setup(str(NML_CHI2K))
    return prob


def test_import_core():
    assert hasattr(RegcoilProblem, "setup")
    assert hasattr(RegcoilProblem, "solve_lambda")


def test_one_lambda_axisymmetry(chi2k_problem):
    chi2_B, chi2_K, max_Bnormal, max_K = chi2k_problem.solve_lambda(1e-15)
    assert abs(chi2_B) < ABS_TOL
    assert abs(max_Bnormal) < ABS_TOL
    assert chi2_K >= 0.0
    assert max_K >= 0.0


def test_solve_ilambda_axisymmetry(chi2k_problem):
    chi2_B, _chi2_K, max_Bnormal, _max_K = chi2k_problem.solve_ilambda(0)
    assert abs(chi2_B) < ABS_TOL
    assert abs(max_Bnormal) < ABS_TOL


def test_two_instances_noninterfering():
    """Two problems with different regularization options must not clobber each other."""
    if not NML_CHI2K.is_file() or not NML_LB.is_file():
        pytest.skip("missing axisymmetry example namelists")

    a = RegcoilProblem()
    b = RegcoilProblem()
    a.set_verbose(False)
    b.set_verbose(False)
    a.setup(str(NML_CHI2K))
    b.setup(str(NML_LB))

    chi2_B_a, chi2_K_a, max_Bn_a, max_K_a = a.solve_lambda(1e-15)
    chi2_B_b, chi2_K_b, max_Bn_b, max_K_b = b.solve_lambda(1e-15)

    assert abs(chi2_B_a) < ABS_TOL
    assert abs(max_Bn_a) < ABS_TOL
    assert abs(chi2_B_b) < ABS_TOL
    assert abs(max_Bn_b) < ABS_TOL

    # Re-solve A after B; results must still be consistent (no shared Fortran globals).
    chi2_B_a2, chi2_K_a2, max_Bn_a2, max_K_a2 = a.solve_lambda(1e-15)
    assert abs(chi2_B_a2 - chi2_B_a) < ABS_TOL
    assert abs(chi2_K_a2 - chi2_K_a) < 1e-8 * max(1.0, abs(chi2_K_a))
    assert abs(max_Bn_a2 - max_Bn_a) < ABS_TOL
    assert abs(max_K_a2 - max_K_a) < 1e-8 * max(1.0, abs(max_K_a))

    assert a.nlambda() == 3
    assert b.nlambda() == 3
    assert chi2_K_a >= 0.0 and chi2_K_b >= 0.0
    assert max_K_a >= 0.0 and max_K_b >= 0.0
