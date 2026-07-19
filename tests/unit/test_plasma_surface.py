"""Unit tests for `PlasmaSurface`.

`from_vmec` is checked against golden r/normal/area/G/curpol values computed
by the legacy compiled Fortran (`regcoil_init_plasma`, geometry_option_plasma=2)
on the same VMEC file -- see tests/unit/_golden.py.
"""

from pathlib import Path

import numpy as np
import pytest

from regcoil import PlasmaSurface

from ._golden import NORMAL_PLASMA, R_PLASMA

REPO_ROOT = Path(__file__).resolve().parents[2]
EQUILIBRIA = REPO_ROOT / "equilibria"
DATA = Path(__file__).resolve().parent / "data"

REF_AREA = 2.4910821007084863e01
REF_VOLUME = 1.6304247536071577e00
REF_G = 1.1884578094260072e07
REF_CURPOL = 4.9782004309255496e00


def test_from_vmec_full_mesh_matches_legacy():
    plasma = PlasmaSurface.from_vmec(str(EQUILIBRIA / "wout_li383_1.4m.nc"), ntheta=4, nzeta=3, mesh="full")

    assert plasma.nfp == 3
    assert plasma.stellarator_symmetric
    np.testing.assert_allclose(plasma.area, REF_AREA, rtol=1e-10)
    np.testing.assert_allclose(plasma.volume, REF_VOLUME, rtol=1e-10)
    np.testing.assert_allclose(plasma.net_poloidal_current_Amperes, REF_G, rtol=1e-10)
    np.testing.assert_allclose(plasma.curpol, REF_CURPOL, rtol=1e-10)
    np.testing.assert_allclose(plasma.r, R_PLASMA, atol=1e-10)
    np.testing.assert_allclose(plasma.normal, NORMAL_PLASMA, atol=1e-10)


def test_from_vmec_mesh_options():
    full = PlasmaSurface.from_vmec(str(EQUILIBRIA / "wout_li383_1.4m.nc"), ntheta=8, nzeta=6, mesh="full")
    half = PlasmaSurface.from_vmec(str(EQUILIBRIA / "wout_li383_1.4m.nc"), ntheta=8, nzeta=6, mesh="half")
    # Different (but nearby) surfaces; both should be well-formed closed tori.
    assert full.area > 0 and half.area > 0
    assert not np.allclose(full.rmnc, half.rmnc)

    with pytest.raises(ValueError):
        PlasmaSurface.from_vmec(str(EQUILIBRIA / "wout_li383_1.4m.nc"), mesh="bogus")


def test_from_vmec_straight_field_line_not_implemented():
    with pytest.raises(NotImplementedError):
        PlasmaSurface.from_vmec(str(EQUILIBRIA / "wout_li383_1.4m.nc"), straight_field_line=True)


def test_from_focus_matches_legacy():
    plasma = PlasmaSurface.from_focus(str(DATA / "focus_boundary_small.txt"), ntheta=6, nzeta=4)

    assert plasma.nfp == 4
    assert plasma.mnmax == 3
    np.testing.assert_allclose(plasma.xm, [0, 1, 1])
    np.testing.assert_allclose(plasma.xn, [0, 0, 4])  # xn scaled by nfp
    np.testing.assert_allclose(plasma.rmnc, [3.0, 0.0, -0.03])
    np.testing.assert_allclose(plasma.rmns, [0.0, 0.0, 0.01])
    np.testing.assert_allclose(plasma.zmnc, [0.0, 0.0, 0.02])
    np.testing.assert_allclose(plasma.zmns, [0.0, 0.5, 0.04])

    np.testing.assert_allclose(plasma.area, 3.9595143548546638e01, rtol=1e-10)
    np.testing.assert_allclose(plasma.volume, 6.8561615968052089e-02, rtol=1e-10)

    bnormal_ref = np.array([
        [1.5000000000000000e-03, -1.5000000000000000e-03, -5.0000000000000023e-04, 2.5000000000000001e-03],
        [1.4378221735089299e-04, -2.8562177826491069e-03, -1.8562177826491073e-03, 1.1437822173508932e-03],
        [1.3562177826491066e-03, -1.6437822173508935e-03, -6.4378221735089365e-04, 2.3562177826491068e-03],
        [1.5000000000000002e-03, -1.4999999999999998e-03, -5.0000000000000001e-04, 2.5000000000000001e-03],
        [1.4378221735089310e-04, -2.8562177826491073e-03, -1.8562177826491073e-03, 1.1437822173508932e-03],
        [1.3562177826491073e-03, -1.6437822173508928e-03, -6.4378221735089311e-04, 2.3562177826491073e-03],
    ])
    np.testing.assert_allclose(plasma.Bnormal_from_plasma_current, bnormal_ref, atol=1e-15)


def test_set_bnormal_from_bnorm_file_matches_legacy():
    plasma = PlasmaSurface.from_vmec(str(EQUILIBRIA / "wout_li383_1.4m.nc"), ntheta=5, nzeta=4, mesh="full")
    np.testing.assert_allclose(plasma.curpol, REF_CURPOL, rtol=1e-10)

    plasma.set_bnormal_from_bnorm_file(str(EQUILIBRIA / "bnorm.li383_1.4m"))

    bnormal_ref = np.array([
        [0.0000000000000000e00, -1.7244761670044260e-02, -2.4703588793970784e-17, 1.7244761670044254e-02],
        [7.3930007084545266e-02, 6.8931311087464781e-02, 6.1293574619278389e-02, 2.3962156630938525e-02],
        [-6.1218251646474209e-02, -6.5383410939484118e-02, 1.0574768071596470e-02, 3.1912723992300787e-02],
        [6.1218251646474202e-02, -3.1912723992300801e-02, -1.0574768071596385e-02, 6.5383410939484105e-02],
        [-7.3930007084545293e-02, -2.3962156630938456e-02, -6.1293574619278396e-02, -6.8931311087464781e-02],
    ])
    np.testing.assert_allclose(plasma.Bnormal_from_plasma_current, bnormal_ref, atol=1e-15)


def test_from_ascii_table_roundtrip(tmp_path):
    # Legacy geometry_option_plasma=6 format: header, mnmax, header, then rows of
    # "m n rmnc zmns rmns zmnc".
    table = tmp_path / "shape.txt"
    table.write_text("comment\n2\ncomment\n0 0 3.0 0.0 0.0 0.0\n1 0 1.5 1.5 0.0 0.0\n")

    plasma = PlasmaSurface.from_ascii_table(str(table), nfp=5, ntheta=8, nzeta=8)

    assert plasma.nfp == 5
    assert not plasma.stellarator_symmetric  # legacy always sets lasym=True for this format
    np.testing.assert_allclose(plasma.xm, [0, 1])
    np.testing.assert_allclose(plasma.xn, [0, 0])
    np.testing.assert_allclose(plasma.rmnc, [3.0, 1.5])
    np.testing.assert_allclose(plasma.zmns, [0.0, 1.5])
    # A plain circular cross-section: area should match circular_torus.
    from regcoil import FourierSurface

    torus = FourierSurface.circular_torus(R0=3.0, a=1.5, nfp=5, ntheta=8, nzeta=8)
    np.testing.assert_allclose(plasma.area, torus.area, rtol=1e-12)
