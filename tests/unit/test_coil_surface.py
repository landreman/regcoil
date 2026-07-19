"""Unit tests for `CoilSurface`.

`from_nescin` is checked against golden r/normal/area/volume values computed
by the legacy compiled Fortran (`regcoil_init_coil_surface`,
geometry_option_coil=3) on the same nescin file -- see tests/unit/_golden.py.
"""

from pathlib import Path

import numpy as np
import pytest

from regcoil import CoilSurface

from ._golden import NORMAL_COIL, R_COIL

REPO_ROOT = Path(__file__).resolve().parents[2]
EQUILIBRIA = REPO_ROOT / "equilibria"

NESCIN_FILE = EQUILIBRIA / "nescin.w7x_standardConfig_separation0.3"

REF_AREA = 1.9033611117049577e02
REF_VOLUME = 4.9883351694634307e01


def test_from_nescin_matches_legacy():
    coil = CoilSurface.from_nescin(str(NESCIN_FILE), nfp=5, ntheta=4, nzeta=3)

    assert coil.mnmax == 163
    np.testing.assert_allclose(coil.area, REF_AREA, rtol=1e-10)
    np.testing.assert_allclose(coil.volume, REF_VOLUME, rtol=1e-10)
    np.testing.assert_allclose(coil.r, R_COIL, atol=1e-9)
    np.testing.assert_allclose(coil.normal, NORMAL_COIL, atol=1e-9)


def test_filter_modes_zeroes_high_harmonics():
    coil = CoilSurface.from_nescin(str(NESCIN_FILE), nfp=5, ntheta=8, nzeta=8)
    assert np.any(np.abs(coil.xm) > 2) or np.any(np.abs(coil.xn) > 2 * coil.nfp)

    coil.filter_modes(mpol_filter=2, ntor_filter=2)
    kept = (np.abs(coil.xm) <= 2) & (np.abs(coil.xn) <= 2 * coil.nfp)
    assert np.all(coil.rmnc[~kept] == 0)
    assert np.all(coil.rmns[~kept] == 0)
    assert np.all(coil.zmnc[~kept] == 0)
    assert np.all(coil.zmns[~kept] == 0)


def test_from_nescin_filter_at_construction():
    filtered = CoilSurface.from_nescin(str(NESCIN_FILE), nfp=5, ntheta=8, nzeta=8, mpol_filter=2, ntor_filter=2)
    unfiltered = CoilSurface.from_nescin(str(NESCIN_FILE), nfp=5, ntheta=8, nzeta=8)
    assert not np.allclose(filtered.rmnc, unfiltered.rmnc)


def test_from_uniform_offset_deferred_to_phase7():
    with pytest.raises(NotImplementedError):
        CoilSurface.from_uniform_offset(plasma=None, separation=0.5)
