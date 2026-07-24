"""Unit tests for VMEC wout I/O, including the optional simsopt `Vmec` path.

simsopt is an optional dependency. Tests that need it are skipped when it
cannot be imported; CI installs simsopt so those tests run there.
"""

from pathlib import Path

import numpy as np
import pytest

import regcoil
from regcoil._io import read_vmec_wout

REPO_ROOT = Path(__file__).resolve().parents[2]
WOUT = REPO_ROOT / "equilibria" / "wout_li383_1.4m.nc"


def test_read_vmec_wout_from_file():
    data = read_vmec_wout(WOUT)
    assert data["nfp"] == 3
    assert data["lasym"] is False
    assert data["rmnc"].ndim == 2
    assert data["xm"].shape == data["xn"].shape


def test_read_vmec_wout_simsopt_vmec_matches_file():
    pytest.importorskip("simsopt")
    from simsopt.mhd import Vmec

    from_file = read_vmec_wout(WOUT)
    from_vmec = read_vmec_wout(Vmec(str(WOUT)))

    assert from_file.keys() == from_vmec.keys()
    for key, expected in from_file.items():
        got = from_vmec[key]
        if expected is None:
            assert got is None
        elif isinstance(expected, (bool, int, np.integer)):
            assert got == expected
        else:
            np.testing.assert_array_equal(got, expected)


def _quickstart_solve(wout, *, ntheta, nzeta, mpol, ntor, lam):
    """Lower-resolution version of the Quickstart example in docs/quickstart.md."""
    ds = regcoil.examples("W7-X")
    plasma = regcoil.PlasmaSurface.from_wout(wout, ntheta=ntheta, nzeta=nzeta)
    plasma.set_bnormal_from_virtual_casing(ds.vcasing)
    coil = regcoil.CoilSurface.from_uniform_offset(
        plasma, separation=0.3, ntheta=ntheta, nzeta=nzeta, mpol=mpol, ntor=ntor,
    )
    problem = regcoil.Regcoil(plasma, coil, mpol_potential=mpol, ntor_potential=ntor)
    return problem.solve(lam=lam)


def test_quickstart_from_simsopt_vmec_matches_wout_file():
    """Quickstart end-to-end: wout path and simsopt `Vmec` must give the same solution."""
    pytest.importorskip("simsopt")
    from simsopt.mhd import Vmec

    ds = regcoil.examples("W7-X")
    kwargs = dict(ntheta=16, nzeta=16, mpol=4, ntor=4, lam=1e-14)

    sol_file = _quickstart_solve(ds.wout, **kwargs)
    sol_vmec = _quickstart_solve(Vmec(str(ds.wout)), **kwargs)

    np.testing.assert_array_equal(sol_file.solution, sol_vmec.solution)
    assert sol_file.f_B == sol_vmec.f_B
    assert sol_file.f_K == sol_vmec.f_K
    assert sol_file.max_K == sol_vmec.max_K
    assert sol_file.max_Bnormal == sol_vmec.max_Bnormal
    assert sol_file.max_Bnormal_over_B == sol_vmec.max_Bnormal_over_B
    assert sol_file.avg_Bnormal_over_B == sol_vmec.avg_Bnormal_over_B
