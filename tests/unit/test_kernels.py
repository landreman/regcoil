"""Unit tests for the Phase 7 stateless Fortran kernels
(`regcoil._core.build_g_and_h`, `.build_inductance`, `.uniform_offset_surface`).

See `_golden_kernels.py` for how the golden reference values were generated
(the legacy compiled Fortran, run standalone outside the package build --
same pattern as `_golden.py` for the Phase 6 surface classes).
"""

import numpy as np
import pytest

import regcoil._core as _core
from regcoil import CoilSurface, PlasmaSurface

from ._golden_kernels import G_GOLDEN, H_GOLDEN, RMNC_OUT, RMNS_OUT, XM_OUT, XN_OUT, ZMNC_OUT, ZMNS_OUT

NT_P, NZ_P = 3, 4
NT_C, NZ_C, NFP = 5, 6, 2
NET_POLOIDAL_CURRENT = 1.0
NET_TOROIDAL_CURRENT = 0.3


def _manufactured_geometry():
    """The deterministic (non-physical) small grid golden.py's g/h case was
    generated from: sin/cos of a linear combination of (component, itheta,
    izeta), well-separated in space (plasma near r~1, coil near r~5) so the
    pairwise kernel never divides by zero.
    """
    r_plasma = np.zeros((3, NT_P, NZ_P), order="F")
    normal_plasma = np.zeros((3, NT_P, NZ_P), order="F")
    for iz in range(1, NZ_P + 1):
        for it in range(1, NT_P + 1):
            for cc in range(1, 4):
                r_plasma[cc - 1, it - 1, iz - 1] = 1.0 + 0.3 * np.sin(0.7 * cc + 0.4 * it + 0.9 * iz)
                normal_plasma[cc - 1, it - 1, iz - 1] = 0.5 + 0.2 * np.cos(0.3 * cc - 0.6 * it + 0.5 * iz)

    nzl = NZ_C * NFP
    r_coil = np.zeros((3, NT_C, nzl), order="F")
    normal_coil = np.zeros((3, NT_C, nzl), order="F")
    drdtheta_coil = np.zeros((3, NT_C, nzl), order="F")
    drdzeta_coil = np.zeros((3, NT_C, nzl), order="F")
    for izl in range(1, nzl + 1):
        for it in range(1, NT_C + 1):
            for cc in range(1, 4):
                r_coil[cc - 1, it - 1, izl - 1] = 5.0 + 0.5 * np.sin(0.5 * cc + 0.3 * it - 0.2 * izl)
                normal_coil[cc - 1, it - 1, izl - 1] = 0.6 + 0.3 * np.cos(0.4 * cc + 0.6 * it + 0.3 * izl)
                drdtheta_coil[cc - 1, it - 1, izl - 1] = 0.3 + 0.15 * np.sin(0.6 * cc - 0.2 * it + 0.4 * izl)
                drdzeta_coil[cc - 1, it - 1, izl - 1] = 0.4 + 0.2 * np.cos(0.2 * cc + 0.5 * it - 0.3 * izl)

    return r_plasma, normal_plasma, r_coil, normal_coil, drdtheta_coil, drdzeta_coil


def _manufactured_basis_functions():
    """mpol_potential=1, ntor_potential=0, symmetry_option=1 (sin only): the
    single basis function is sin(theta_coil), matching the legacy
    `regcoil_build_matrices` basis-function convention.
    """
    theta_coil = 2 * np.pi * np.arange(NT_C) / NT_C
    basis_functions = np.zeros((NT_C * NZ_C, 1), order="F")
    for izeta_coil in range(NZ_C):
        for itheta_coil in range(NT_C):
            index = izeta_coil * NT_C + itheta_coil
            basis_functions[index, 0] = np.sin(theta_coil[itheta_coil])
    return basis_functions


def _dtheta_dzeta_coil():
    return 2 * np.pi / NT_C, (2 * np.pi / NFP) / NZ_C


def test_build_g_and_h_matches_build_inductance_times_basis():
    r_plasma, normal_plasma, r_coil, normal_coil, drdtheta_coil, drdzeta_coil = _manufactured_geometry()
    basis_functions = _manufactured_basis_functions()
    dtheta_coil, dzeta_coil = _dtheta_dzeta_coil()

    g, h = _core.build_g_and_h(
        r_plasma, normal_plasma, r_coil, normal_coil, drdtheta_coil, drdzeta_coil,
        basis_functions, NFP, NET_POLOIDAL_CURRENT, NET_TOROIDAL_CURRENT, dtheta_coil, dzeta_coil,
    )
    inductance, h_from_inductance = _core.build_inductance(
        r_plasma, normal_plasma, r_coil, normal_coil, drdtheta_coil, drdzeta_coil,
        NFP, NET_POLOIDAL_CURRENT, NET_TOROIDAL_CURRENT, dtheta_coil, dzeta_coil,
    )
    g_direct = dtheta_coil * dzeta_coil * (inductance @ basis_functions)

    np.testing.assert_allclose(g, g_direct, atol=1e-18, rtol=1e-12)
    np.testing.assert_allclose(h, h_from_inductance, atol=1e-24, rtol=1e-12)


def test_build_g_and_h_matches_legacy_golden():
    r_plasma, normal_plasma, r_coil, normal_coil, drdtheta_coil, drdzeta_coil = _manufactured_geometry()
    basis_functions = _manufactured_basis_functions()
    dtheta_coil, dzeta_coil = _dtheta_dzeta_coil()

    g, h = _core.build_g_and_h(
        r_plasma, normal_plasma, r_coil, normal_coil, drdtheta_coil, drdzeta_coil,
        basis_functions, NFP, NET_POLOIDAL_CURRENT, NET_TOROIDAL_CURRENT, dtheta_coil, dzeta_coil,
    )

    np.testing.assert_allclose(g.ravel(order="F"), G_GOLDEN, rtol=1e-10)
    np.testing.assert_allclose(h, H_GOLDEN, rtol=1e-10)


def test_build_inductance_matches_legacy_golden():
    r_plasma, normal_plasma, r_coil, normal_coil, drdtheta_coil, drdzeta_coil = _manufactured_geometry()
    basis_functions = _manufactured_basis_functions()
    dtheta_coil, dzeta_coil = _dtheta_dzeta_coil()

    inductance, h = _core.build_inductance(
        r_plasma, normal_plasma, r_coil, normal_coil, drdtheta_coil, drdzeta_coil,
        NFP, NET_POLOIDAL_CURRENT, NET_TOROIDAL_CURRENT, dtheta_coil, dzeta_coil,
    )
    g_direct = dtheta_coil * dzeta_coil * (inductance @ basis_functions)

    np.testing.assert_allclose(g_direct.ravel(order="F"), G_GOLDEN, rtol=1e-10)
    np.testing.assert_allclose(h, H_GOLDEN, rtol=1e-10)


def test_kernels_reject_non_fortran_contiguous_arrays():
    r_plasma, normal_plasma, r_coil, normal_coil, drdtheta_coil, drdzeta_coil = _manufactured_geometry()
    basis_functions = _manufactured_basis_functions()
    dtheta_coil, dzeta_coil = _dtheta_dzeta_coil()

    with pytest.raises(ValueError, match="Fortran-contiguous"):
        _core.build_g_and_h(
            np.ascontiguousarray(r_plasma), normal_plasma, r_coil, normal_coil, drdtheta_coil, drdzeta_coil,
            basis_functions, NFP, NET_POLOIDAL_CURRENT, NET_TOROIDAL_CURRENT, dtheta_coil, dzeta_coil,
        )


def test_kernels_reject_wrong_dtype():
    r_plasma, normal_plasma, r_coil, normal_coil, drdtheta_coil, drdzeta_coil = _manufactured_geometry()
    dtheta_coil, dzeta_coil = _dtheta_dzeta_coil()

    with pytest.raises(TypeError, match="float64"):
        _core.build_inductance(
            r_plasma.astype(np.float32), normal_plasma, r_coil, normal_coil, drdtheta_coil, drdzeta_coil,
            NFP, NET_POLOIDAL_CURRENT, NET_TOROIDAL_CURRENT, dtheta_coil, dzeta_coil,
        )


def test_two_instances_different_sizes_do_not_interfere():
    """Stateless kernel: repeated/interleaved calls at different problem
    sizes must not corrupt each other (no persistent Fortran state)."""
    r_plasma, normal_plasma, r_coil, normal_coil, drdtheta_coil, drdzeta_coil = _manufactured_geometry()
    basis_functions = _manufactured_basis_functions()
    dtheta_coil, dzeta_coil = _dtheta_dzeta_coil()

    rng = np.random.default_rng(1)
    r_plasma_b = np.asfortranarray(rng.normal(size=(3, 2, 3)) * 0.2 + 1.5)
    normal_plasma_b = np.asfortranarray(rng.normal(size=(3, 2, 3)))
    r_coil_b = np.asfortranarray(rng.normal(size=(3, 3, 4 * 3)) * 0.3 + 6.0)
    normal_coil_b = np.asfortranarray(rng.normal(size=(3, 3, 4 * 3)))
    drdtheta_coil_b = np.asfortranarray(rng.normal(size=(3, 3, 4 * 3)))
    drdzeta_coil_b = np.asfortranarray(rng.normal(size=(3, 3, 4 * 3)))
    basis_functions_b = np.asfortranarray(rng.normal(size=(3 * 4, 2)))

    g_a1, h_a1 = _core.build_g_and_h(
        r_plasma, normal_plasma, r_coil, normal_coil, drdtheta_coil, drdzeta_coil,
        basis_functions, NFP, NET_POLOIDAL_CURRENT, NET_TOROIDAL_CURRENT, dtheta_coil, dzeta_coil,
    )
    g_b, h_b = _core.build_g_and_h(
        r_plasma_b, normal_plasma_b, r_coil_b, normal_coil_b, drdtheta_coil_b, drdzeta_coil_b,
        basis_functions_b, 3, 0.5, 0.1, 0.2, 0.3,
    )
    g_a2, h_a2 = _core.build_g_and_h(
        r_plasma, normal_plasma, r_coil, normal_coil, drdtheta_coil, drdzeta_coil,
        basis_functions, NFP, NET_POLOIDAL_CURRENT, NET_TOROIDAL_CURRENT, dtheta_coil, dzeta_coil,
    )

    np.testing.assert_array_equal(g_a1, g_a2)
    np.testing.assert_array_equal(h_a1, h_a2)
    assert g_b.shape == (2 * 3, 2)
    assert not np.allclose(g_b.shape, g_a1.shape)


def test_uniform_offset_surface_matches_legacy_golden():
    """Plasma with an (m=1, n=1) helical bump, so the offset surface is not a
    circular torus: this forces a genuine toroidal-angle root solve
    (`atan2(y, x)` is not linear in the root-solve variable, unlike the
    circular-cross-section case where the residual is identically zero).

    The golden values come from the legacy compiled Fortran. Since ADR-031 the
    root solve is numpy (`_solve_zeta_for_phi`, Illinois-modified regula falsi)
    rather than `regcoil_fzero`, so this is the direct legacy check on the
    ported implementation.
    """
    plasma = PlasmaSurface(
        xm=[0, 1, 1], xn=[0, 0, 3],  # xn already includes the nfp=3 factor
        rmnc=[5.0, 1.0, 0.15], zmns=[0.0, 1.0, 0.1],
        nfp=3, ntheta=8, nzeta=8,
    )
    coil = CoilSurface.from_uniform_offset(
        plasma, separation=0.4, ntheta=8, nzeta=8, mpol=3, ntor=2,
        standard_toroidal_angle=True, ntheta_transform=6, nzeta_transform=5, tol=1e-10,
    )

    np.testing.assert_array_equal(coil.xm, XM_OUT)
    np.testing.assert_array_equal(coil.xn, XN_OUT)
    np.testing.assert_allclose(coil.rmnc, RMNC_OUT, atol=1e-12, rtol=1e-10)
    np.testing.assert_allclose(coil.rmns, RMNS_OUT, atol=1e-12)
    np.testing.assert_allclose(coil.zmnc, ZMNC_OUT, atol=1e-12)
    np.testing.assert_allclose(coil.zmns, ZMNS_OUT, atol=1e-12, rtol=1e-10)


def test_uniform_offset_surface_circular_torus_is_exact():
    """A circular-cross-section torus offset outward by `separation` stays a
    circular-cross-section torus of minor radius a + separation -- an exact
    analytic check independent of the golden legacy values above."""
    plasma = PlasmaSurface.circular_torus(R0=5.0, a=1.0, nfp=3, ntheta=8, nzeta=8)
    coil = CoilSurface.from_uniform_offset(
        plasma, separation=0.3, ntheta=8, nzeta=8, mpol=2, ntor=1,
        standard_toroidal_angle=True, ntheta_transform=16, nzeta_transform=16,
    )

    for m, n, rc, zs in zip(coil.xm, coil.xn, coil.rmnc, coil.zmns):
        if m == 0 and n == 0:
            assert rc == pytest.approx(5.0, abs=1e-9)
        elif m == 1 and n == 0:
            assert rc == pytest.approx(1.3, abs=1e-9)
            assert zs == pytest.approx(1.3, abs=1e-9)
        else:
            assert rc == pytest.approx(0.0, abs=1e-8)
            assert zs == pytest.approx(0.0, abs=1e-8)
