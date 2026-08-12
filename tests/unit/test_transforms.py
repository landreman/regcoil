"""Unit tests for `regcoil._transforms`, the separated double-Fourier sums.

Each routine there replaces a direct sum over modes with one that factors the
angle `m*theta - n*zeta` into its two halves. The factorization is exact, so
every test here is the same statement: the fast form equals the direct one.
"""

import numpy as np
import pytest

from regcoil._transforms import (
    dft_weight,
    evaluate_modes,
    evaluate_modes_by_column,
    fit_modes,
)
from regcoil.coil_surface import _uniform_offset_modes

NFP = 3


def _modes_and_amplitudes(mpol=4, ntor=3, nfields=1, seed=0):
    xm, xn = _uniform_offset_modes(NFP, mpol, ntor)
    rng = np.random.default_rng(seed)
    shape = (nfields, xm.size) if nfields > 1 else (xm.size,)
    return xm, xn, rng.normal(size=shape), rng.normal(size=shape)


def _direct_sum(cmn, smn, xm, xn, theta, zeta):
    """The definition: build every mode's angle and add them up."""
    angle = xm[:, None, None] * theta[None, :, None] - xn[:, None, None] * zeta[None, None, :]
    return np.tensordot(cmn, np.cos(angle), axes=(-1, 0)) + np.tensordot(
        smn, np.sin(angle), axes=(-1, 0)
    )


@pytest.mark.parametrize("nfields", [1, 3])
def test_evaluate_modes_matches_the_direct_sum(nfields):
    xm, xn, cmn, smn = _modes_and_amplitudes(nfields=nfields)
    # Deliberately non-uniform grids: nothing in the factorization needs them
    # to be uniform.
    theta = np.array([0.0, 0.3, 1.1, 2.0, 5.9])
    zeta = np.array([0.0, 0.2, 0.7])

    np.testing.assert_allclose(
        evaluate_modes(cmn, smn, xm, xn, theta, zeta),
        _direct_sum(cmn, smn, xm, xn, theta, zeta),
        atol=1e-13,
    )


def test_evaluate_modes_by_column_matches_the_direct_sum():
    """The column form's theta may differ from one zeta to the next, so check
    it against a point-by-point evaluation of the series."""
    xm, xn, cmn, smn = _modes_and_amplitudes(nfields=3, seed=1)
    zeta = np.array([0.0, 0.2, 0.7, 1.4])
    rng = np.random.default_rng(2)
    theta = np.sort(rng.uniform(0, 2 * np.pi, (6, zeta.size)), axis=0)

    result = evaluate_modes_by_column(cmn, smn, xm, xn, theta, zeta)

    expected = np.empty_like(result)
    for i in range(theta.shape[0]):
        for j in range(zeta.size):
            angle = xm * theta[i, j] - xn * zeta[j]
            expected[:, i, j] = cmn @ np.cos(angle) + smn @ np.sin(angle)
    np.testing.assert_allclose(result, expected, atol=1e-13)


def test_evaluate_modes_by_column_reduces_to_the_grid_form():
    xm, xn, cmn, smn = _modes_and_amplitudes(seed=3)
    theta = np.array([0.0, 0.9, 2.2, 4.5])
    zeta = np.array([0.1, 0.6, 1.0])

    np.testing.assert_allclose(
        evaluate_modes_by_column(
            cmn, smn, xm, xn, np.broadcast_to(theta[:, None], (4, 3)), zeta
        ),
        evaluate_modes(cmn, smn, xm, xn, theta, zeta),
        atol=1e-14,
    )


def test_fit_modes_matches_the_direct_sum():
    xm, xn = _uniform_offset_modes(NFP, 4, 3)
    theta = 2 * np.pi * np.arange(16) / 16
    zeta = (2 * np.pi / NFP) * np.arange(9) / 9
    field = np.random.default_rng(4).normal(size=(16, 9))

    cmn, smn = fit_modes(field, xm, xn, theta, zeta)

    angle = xm[:, None, None] * theta[None, :, None] - xn[:, None, None] * zeta[None, None, :]
    np.testing.assert_allclose(
        cmn, np.tensordot(np.cos(angle), field, axes=([1, 2], [0, 1])), atol=1e-11
    )
    np.testing.assert_allclose(
        smn, np.tensordot(np.sin(angle), field, axes=([1, 2], [0, 1])), atol=1e-11
    )


@pytest.mark.parametrize("ntheta,nzeta", [(16, 8), (15, 9), (16, 9)])
def test_fit_then_evaluate_is_the_identity_at_the_nyquist_limit(ntheta, nzeta):
    """What `dft_weight`'s halved Nyquist weight buys: transforming a field
    sampled on the grid and evaluating the result back on it returns the field,
    when the modes reach the grid's Nyquist limit."""
    xm, xn = _uniform_offset_modes(NFP, ntheta // 2, nzeta // 2)
    theta = 2 * np.pi * np.arange(ntheta) / ntheta
    zeta = (2 * np.pi / NFP) * np.arange(nzeta) / nzeta
    field = np.random.default_rng(5).normal(size=(ntheta, nzeta))

    cmn, smn = fit_modes(field, xm, xn, theta, zeta)
    weight = dft_weight(xm, xn, ntheta, nzeta, NFP)
    cmn, smn = cmn * weight, smn * weight
    cmn[0] = field.mean()  # the m = n = 0 amplitude is the average, not twice it
    smn[0] = 0.0

    np.testing.assert_allclose(
        evaluate_modes(cmn, smn, xm, xn, theta, zeta), field, atol=1e-12
    )
