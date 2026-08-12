"""Double-Fourier series on a `(theta, zeta)` grid, evaluated and fitted
without ever forming a per-mode 3D array.

Every sum here is over modes `(m, n)` in the package's convention (angle
`m*theta - n*zeta`, with `xn` already carrying the factor of `nfp`):

    f(theta, zeta) = sum_mn cmn*cos(m*theta - n*zeta) + smn*sin(m*theta - n*zeta)

The direct way to evaluate that on an `ntheta x nzeta` grid is to build the
`(mnmax, ntheta, nzeta)` array of angles and reduce it. Expanding the angle
difference instead,

    cos(m*theta - n*zeta) = cos(m*theta)*cos(n*zeta) + sin(m*theta)*sin(n*zeta)
    sin(m*theta - n*zeta) = sin(m*theta)*cos(n*zeta) - cos(m*theta)*sin(n*zeta)

separates the two angles, so the sum over `n` can be taken first, leaving
`nm` poloidal amplitudes per `zeta` (`poloidal_amplitudes`) that the sum over
`m` then combines. Here `nm`/`nn` are the numbers of *distinct* `m`/`n` in the
mode list -- for a full `mpol x ntor` list, `nm + nn` values of the trig
functions in place of `mnmax ~ nm*nn` of them, evaluated on the grid.

Two grid shapes are supported, and the split above is what lets the second one
be nearly as cheap as the first:

* `evaluate_modes`, the tensor-product grid `theta x zeta`, where the sum over
  `m` is a matrix product;
* `evaluate_modes_by_column`, `theta` varying from one `zeta` to the next --
  the shape a theta reparameterization produces -- where it is a contraction
  over the (short) `m` axis of an `(nm, ntheta, nzeta)` array. Reaching for
  the general `Surface.evaluate_at` there instead would cost a factor of `nn`
  in both trig evaluations and memory.

Derivatives need no separate machinery: differentiating the series above shows
that `df/dtheta` is the same series with amplitudes `(cmn, smn) -> (m*smn,
-m*cmn)`, and `df/dzeta` with `(cmn, smn) -> (-n*smn, n*cmn)`, so a caller
that wants a field and its derivatives passes them as extra rows of one
stacked amplitude array (see `FourierSurface._evaluate`).

Nothing here assumes the grids are uniform, or that the mode list fills the
`nm x nn` rectangle: the rectangle is padded with zeros where a mode is
missing. A mode list far sparser than a rectangle (`mnmax << nm*nn`) would
make the padding wasteful rather than wrong; the lists this package sees
(`_uniform_offset_modes`, VMEC's, and nescin's) are rectangular by
construction.
"""

from __future__ import annotations

import numpy as np


def _rectangle(cmn, smn, xm, xn):
    """Scatter the amplitudes onto the `(nm, nn)` rectangle of distinct
    `(m, n)` values, so that the sums over `m` and over `n` can be taken
    independently.

    `cmn`/`smn` are `(mnmax,)` or `(nfields, mnmax)`. Returns
    `(m_values, n_values, rect_c, rect_s)` with the rectangles shaped
    `(nfields, nm, nn)`.
    """
    cmn = np.atleast_2d(cmn)
    smn = np.atleast_2d(smn)
    m_values, m_index = np.unique(xm, return_inverse=True)
    n_values, n_index = np.unique(xn, return_inverse=True)
    m_index = np.ravel(m_index)
    n_index = np.ravel(n_index)

    shape = (cmn.shape[0], m_values.size, n_values.size)
    rect_c = np.zeros(shape)
    rect_s = np.zeros(shape)
    # `add.at` rather than plain assignment so that a mode list which happens
    # to repeat an (m, n) pair still sums, as the series says it should.
    np.add.at(rect_c, (slice(None), m_index, n_index), cmn)
    np.add.at(rect_s, (slice(None), m_index, n_index), smn)
    return m_values.astype(np.float64), n_values.astype(np.float64), rect_c, rect_s


def poloidal_amplitudes(cmn, smn, xm, xn, zeta):
    """Take the sum over `n`, leaving a series in `theta` alone at each `zeta`:

        f(theta, zeta_j) = sum_a amp_c[a, j]*cos(m_a*theta)
                         + amp_s[a, j]*sin(m_a*theta)

    Returns `(m_values, amp_c, amp_s)`, the amplitudes shaped
    `(nfields, nm, len(zeta))`.
    """
    m_values, n_values, rect_c, rect_s = _rectangle(cmn, smn, xm, xn)
    n_zeta = np.outer(n_values, zeta)
    cos_n = np.cos(n_zeta)
    sin_n = np.sin(n_zeta, out=n_zeta)  # `n_zeta` is dead after this
    return m_values, rect_c @ cos_n - rect_s @ sin_n, rect_c @ sin_n + rect_s @ cos_n


def evaluate_modes(cmn, smn, xm, xn, theta, zeta):
    """Sum the series on the tensor-product grid `theta x zeta`.

    `cmn`/`smn` are the cosine/sine amplitudes of the modes `(xm, xn)`, either
    `(mnmax,)` for a single field or `(nfields, mnmax)` for several sharing the
    mode list. The result is `(len(theta), len(zeta))`, or
    `(nfields, len(theta), len(zeta))` in the stacked case.
    """
    m_values, amp_c, amp_s = poloidal_amplitudes(cmn, smn, xm, xn, zeta)
    m_theta = np.outer(m_values, theta)
    cos_m = np.cos(m_theta)
    sin_m = np.sin(m_theta, out=m_theta)
    values = cos_m.T @ amp_c + sin_m.T @ amp_s
    return values if np.ndim(cmn) > 1 else values[0]


def evaluate_modes_by_column(cmn, smn, xm, xn, theta, zeta):
    """Sum the series at `(theta[i, j], zeta[j])`: one `zeta` grid, but a
    poloidal angle that may differ from column to column.

    `theta` is `(ntheta, len(zeta))`. Amplitudes and result are shaped as in
    `evaluate_modes`.
    """
    theta = np.asarray(theta, dtype=np.float64)
    m_values, amp_c, amp_s = poloidal_amplitudes(cmn, smn, xm, xn, zeta)
    m_theta = m_values[:, None, None] * theta[None, :, :]
    cos_m = np.cos(m_theta)
    sin_m = np.sin(m_theta, out=m_theta)
    values = np.einsum("fmj,mij->fij", amp_c, cos_m, optimize=True)
    values += np.einsum("fmj,mij->fij", amp_s, sin_m, optimize=True)
    return values if np.ndim(cmn) > 1 else values[0]


def fit_modes(field, xm, xn, theta, zeta):
    """The unweighted discrete-transform sums of `field` against each mode,

        cmn_k = sum_ij cos(xm_k*theta_i - xn_k*zeta_j) * field_ij,

    and likewise `smn_k` with the sine -- the transform `evaluate_modes`
    inverts. `field` is `(len(theta), len(zeta))` and the returned arrays are
    `(mnmax,)`; the caller supplies the normalization (see `dft_weight`).
    """
    m_values, m_index = np.unique(xm, return_inverse=True)
    n_values, n_index = np.unique(xn, return_inverse=True)
    m_theta = np.outer(m_values.astype(np.float64), theta)
    n_zeta = np.outer(n_values.astype(np.float64), zeta)

    # Contract over theta first, then over zeta, so that the (mnmax, ntheta,
    # nzeta) array of the direct sum never exists.
    field = np.ascontiguousarray(field)
    partial_c = np.cos(m_theta) @ field
    partial_s = np.sin(m_theta, out=m_theta) @ field
    cos_n = np.cos(n_zeta)
    sin_n = np.sin(n_zeta, out=n_zeta)
    rect_c = partial_c @ cos_n.T + partial_s @ sin_n.T
    rect_s = partial_s @ cos_n.T - partial_c @ sin_n.T
    return rect_c[np.ravel(m_index), np.ravel(n_index)], rect_s[np.ravel(m_index), np.ravel(n_index)]


def dft_weight(xm, xn, ntheta, nzeta, nfp):
    """Normalization turning `fit_modes`' raw sums into Fourier amplitudes:
    `2 / (ntheta*nzeta)`, halved for a Nyquist mode in either angle so that
    evaluating the fitted series back on the same grid is the identity.

    The `m = n = 0` amplitude is *not* handled here: it is the plain average
    rather than twice it, and callers set it explicitly.
    """
    weight = np.full(np.shape(xm), 2.0 / (ntheta * nzeta))
    if ntheta % 2 == 0:
        weight[xm == ntheta // 2] /= 2
    if nzeta % 2 == 0:
        weight[np.abs(xn) == nfp * (nzeta // 2)] /= 2
    return weight
