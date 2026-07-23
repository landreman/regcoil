"""File readers for REGCOIL geometry/physics input formats.

Pure Python (no Fortran). VMEC `wout` reading tries
`scipy.io.netcdf_file` first (already a required dependency, and sufficient
for classic-format wout files) and falls back to the optional `netCDF4`
package for HDF5-based files.
"""

from __future__ import annotations

import os

import numpy as np


def _open_netcdf(filename):
    """Open a NetCDF file and return `(variables, get)`.

    `get(name)` returns the named variable as a plain numpy array, for either
    backend and for any rank (including 0-d scalars).
    """
    filename = str(filename)
    try:
        from scipy.io import netcdf_file

        f = netcdf_file(filename, "r", mmap=False)
        variables = f.variables

        def get(name):
            # scipy's netcdf_variable exposes the raw array via `.data` for
            # any rank; `[:]` fails on 0-d (scalar) variables.
            return np.asarray(variables[name].data).copy()

    except Exception as scipy_exc:
        try:
            import netCDF4
        except ImportError:
            raise OSError(
                f"Could not read {filename!r} as a classic NetCDF file with scipy, and "
                "netCDF4 is not installed to try the HDF5-based format. Install the "
                "optional 'netcdf' extra (pip install regcoil[netcdf])."
            ) from scipy_exc
        f = netCDF4.Dataset(filename, "r")
        variables = f.variables

        def get(name):
            return np.asarray(variables[name][:])

    return variables, get


def read_vmec_wout(filename):
    """Read the subset of a VMEC `wout` file that REGCOIL needs."""
    variables, get = _open_netcdf(filename)

    nfp = int(get("nfp"))
    if "lasym__logical__" in variables:
        lasym = bool(get("lasym__logical__"))
    else:
        lasym = bool(get("lasym"))

    xm = get("xm").astype(np.int64)
    xn = get("xn").astype(np.int64)
    # |B| and other Nyquist quantities may use a denser mode set than R/Z.
    if "xm_nyq" in variables:
        xm_nyq = get("xm_nyq").astype(np.int64)
        xn_nyq = get("xn_nyq").astype(np.int64)
    else:
        xm_nyq, xn_nyq = xm, xn
    rmnc = get("rmnc")
    zmns = get("zmns")
    rmns = get("rmns") if lasym and "rmns" in variables else None
    zmnc = get("zmnc") if lasym and "zmnc" in variables else None
    bmnc = get("bmnc")
    bmns = get("bmns") if lasym and "bmns" in variables else None
    bvco = get("bvco")
    bsubvmnc = get("bsubvmnc")
    lmns = get("lmns") if "lmns" in variables else None
    mpol = int(get("mpol"))
    ntor = int(get("ntor"))

    return dict(
        nfp=nfp, lasym=lasym, xm=xm, xn=xn, xm_nyq=xm_nyq, xn_nyq=xn_nyq,
        rmnc=rmnc, zmns=zmns, rmns=rmns, zmnc=zmnc, bmnc=bmnc, bmns=bmns,
        bvco=bvco, bsubvmnc=bsubvmnc, lmns=lmns, mpol=mpol, ntor=ntor,
    )


def read_nescin(filename, nfp):
    """Read a nescin-format coil surface file.

    nescin uses cos(m*theta + n*nfp*zeta); convert xn to the m*theta - n*zeta
    convention used throughout this package.
    """
    match_string = "------ Current Surface"
    with open(filename) as fh:
        lines = fh.readlines()

    for i, line in enumerate(lines):
        if line.startswith(match_string):
            break
    else:
        raise ValueError(f"Could not find {match_string!r} in {filename}")

    mnmax = int(lines[i + 2].split()[0])
    start = i + 5

    xm = np.zeros(mnmax, dtype=np.int64)
    xn = np.zeros(mnmax, dtype=np.int64)
    rmnc = np.zeros(mnmax)
    zmns = np.zeros(mnmax)
    rmns = np.zeros(mnmax)
    zmnc = np.zeros(mnmax)
    for k in range(mnmax):
        vals = lines[start + k].split()
        xm[k] = int(vals[0])
        xn[k] = int(vals[1])
        rmnc[k], zmns[k], rmns[k], zmnc[k] = (float(v) for v in vals[2:6])

    xn = -nfp * xn
    return dict(xm=xm, xn=xn, rmnc=rmnc, zmns=zmns, rmns=rmns, zmnc=zmnc)


def read_ascii_table(filename):
    """Read the plain ASCII-table plasma-shape format (legacy geometry_option_plasma=6)."""
    with open(filename) as fh:
        lines = fh.readlines()

    mnmax = int(lines[1].split()[0])
    start = 3

    xm = np.zeros(mnmax, dtype=np.int64)
    xn = np.zeros(mnmax, dtype=np.int64)
    rmnc = np.zeros(mnmax)
    zmns = np.zeros(mnmax)
    rmns = np.zeros(mnmax)
    zmnc = np.zeros(mnmax)
    for k in range(mnmax):
        vals = lines[start + k].split()
        xm[k] = int(vals[0])
        xn[k] = int(vals[1])
        rmnc[k], zmns[k], rmns[k], zmnc[k] = (float(v) for v in vals[2:6])

    return dict(xm=xm, xn=xn, rmnc=rmnc, zmns=zmns, rmns=rmns, zmnc=zmnc)


def read_focus_boundary(filename):
    """Read a FOCUS-format plasma boundary file (`n m Rbc Rbs Zbc Zbs` + optional Bn modes).

    See https://princetonuniversity.github.io/FOCUS/rdsurf.pdf. `xn`/`bfn` are
    converted to the `nfp`-included convention used throughout this package.
    """
    with open(filename) as fh:
        lines = fh.readlines()

    mnmax, nfp, nbf = (int(x) for x in lines[1].split()[:3])
    start = 4

    xm = np.zeros(mnmax, dtype=np.int64)
    xn = np.zeros(mnmax, dtype=np.int64)
    rmnc = np.zeros(mnmax)
    rmns = np.zeros(mnmax)
    zmnc = np.zeros(mnmax)
    zmns = np.zeros(mnmax)
    for i in range(mnmax):
        vals = lines[start + i].split()
        xn[i] = int(vals[0])
        xm[i] = int(vals[1])
        rmnc[i], rmns[i], zmnc[i], zmns[i] = (float(v) for v in vals[2:6])
    xn = xn * nfp

    bfm = bfn = bfc = bfs = None
    if nbf > 0:
        bstart = start + mnmax + 2
        bfn = np.zeros(nbf, dtype=np.int64)
        bfm = np.zeros(nbf, dtype=np.int64)
        bfc = np.zeros(nbf)
        bfs = np.zeros(nbf)
        for i in range(nbf):
            vals = lines[bstart + i].split()
            bfn[i] = int(vals[0])
            bfm[i] = int(vals[1])
            bfc[i], bfs[i] = float(vals[2]), float(vals[3])
        bfn = bfn * nfp

    return dict(
        mnmax=mnmax, nfp=nfp, nbf=nbf, xm=xm, xn=xn, rmnc=rmnc, rmns=rmns, zmnc=zmnc, zmns=zmns,
        bfm=bfm, bfn=bfn, bfc=bfc, bfs=bfs,
    )


def read_bnorm_file(filename, theta, zeta, nfp, curpol):
    """Read a BNORM-format file and evaluate B_normal on the given (theta, zeta) grid.

    BNORM scales B_n by curpol = (2*pi/nfp)*bsubv(m=0,n=0) (the extrapolation
    to the last full VMEC mesh point); undo that scaling with `curpol`.
    """
    Bnormal = np.zeros((len(theta), len(zeta)))
    with open(filename) as fh:
        for line in fh:
            parts = line.split()
            if not parts:
                continue
            m, n, bf = int(parts[0]), int(parts[1]), float(parts[2])
            angle = m * theta[:, None] + n * nfp * zeta[None, :]
            Bnormal += bf * np.sin(angle)
    return Bnormal * curpol


def read_virtual_casing_file(filename):
    """Read a simsopt virtual-casing NetCDF file (`vcasing*.nc`).

    Only the fields REGCOIL needs are returned: the target grid and
    `B_external_normal` on it. Reading is done directly with NetCDF, so simsopt
    need not be installed.
    """
    _, get = _open_netcdf(filename)
    return dict(
        nfp=int(get("nfp")),
        trgt_theta=np.asarray(get("trgt_theta"), dtype=float).ravel(),
        trgt_phi=np.asarray(get("trgt_phi"), dtype=float).ravel(),
        B_external_normal=np.asarray(get("B_external_normal"), dtype=float),
    )


def virtual_casing_data(source):
    """Normalize a simsopt virtual-casing result to the dict
    `read_virtual_casing_file` returns.

    `source` is either a path to a `vcasing*.nc` file, or a
    `simsopt.mhd.VirtualCasing`-like object (anything carrying `nfp`,
    `trgt_theta`, `trgt_phi`, and `B_external_normal` attributes -- note that
    `VirtualCasing.load` leaves scalars as 0-d arrays, which is handled here).
    """
    if isinstance(source, (str, os.PathLike)):
        return read_virtual_casing_file(source)

    missing = [
        name
        for name in ("nfp", "trgt_theta", "trgt_phi", "B_external_normal")
        if not hasattr(source, name)
    ]
    if missing:
        raise TypeError(
            f"{type(source).__name__} is neither a path to a simsopt virtual-casing "
            f"NetCDF file nor a VirtualCasing-like object: missing attribute(s) "
            f"{', '.join(missing)}."
        )
    return dict(
        nfp=int(np.asarray(source.nfp).ravel()[0]),
        trgt_theta=np.asarray(source.trgt_theta, dtype=float).ravel(),
        trgt_phi=np.asarray(source.trgt_phi, dtype=float).ravel(),
        B_external_normal=np.asarray(source.B_external_normal, dtype=float),
    )


def _uniform_grid_start(points, period, name):
    """Check that `points` is a uniformly spaced grid covering `[start, start+period)`,
    and return `start`. Raises if the grid is not uniform or does not tile the period.
    """
    n = len(points)
    if n < 2:
        raise ValueError(f"Need at least 2 {name} grid points, got {n}.")
    spacing = period / n
    expected = points[0] + spacing * np.arange(n)
    if not np.allclose(points, expected, rtol=0, atol=1e-9 * period):
        raise ValueError(
            f"The virtual-casing {name} grid is not uniformly spaced over a span of "
            f"{period}, as REGCOIL's spectral interpolation requires."
        )
    return float(points[0])


def _spectral_interpolate(f, u_targets, v_targets):
    """Trigonometric interpolation of `f` off its sampling grid.

    `f[j, i]` is sampled at `u = j / f.shape[0]`, `v = i / f.shape[1]`, with both
    `u` and `v` periodic with period 1. Returns `f` evaluated on the tensor
    product of `u_targets` and `v_targets`, shape `(len(u_targets), len(v_targets))`.

    This is exact for data band-limited below the Nyquist frequency of the
    sampling grid, and reproduces `f` exactly at the sample points. For an even
    number of samples the Nyquist mode is split evenly between `+N/2` and
    `-N/2`, the standard choice that keeps the interpolant real and even about
    the sample points.
    """
    nu, nv = f.shape
    coeffs = np.fft.fft2(f) / (nu * nv)
    ku = np.fft.fftfreq(nu, d=1.0 / nu)
    kv = np.fft.fftfreq(nv, d=1.0 / nv)
    if nu % 2 == 0:
        j = nu // 2
        coeffs = np.concatenate([coeffs, 0.5 * coeffs[[j]]], axis=0)
        coeffs[j] *= 0.5
        ku = np.concatenate([ku, [-ku[j]]])
    if nv % 2 == 0:
        j = nv // 2
        coeffs = np.concatenate([coeffs, 0.5 * coeffs[:, [j]]], axis=1)
        coeffs[:, j] *= 0.5
        kv = np.concatenate([kv, [-kv[j]]])

    basis_u = np.exp(2j * np.pi * np.outer(np.asarray(u_targets, dtype=float), ku))
    basis_v = np.exp(2j * np.pi * np.outer(np.asarray(v_targets, dtype=float), kv))
    return np.real(basis_u @ coeffs @ basis_v.T)


def bnormal_from_virtual_casing(data, theta, zeta, nfp):
    """Evaluate a simsopt virtual-casing `B_external_normal` on a REGCOIL
    `(theta, zeta)` grid.

    `data` is the dict returned by `virtual_casing_data`. `theta` and `zeta` are
    REGCOIL's angle grids in radians (`zeta` covering one field period); simsopt
    uses the same VMEC poloidal angle and the same standard toroidal angle, but
    normalized to a period of 1 rather than 2*pi.

    simsopt supplies `B_external_normal` on its own fixed uniform grid, which in
    general differs from REGCOIL's, so the data are interpolated. The
    interpolation is trigonometric (`_spectral_interpolate`) rather than
    polynomial: the data are periodic and smooth, so this converges spectrally
    and is exact where the two grids happen to coincide.

    Two simsopt grid layouts are supported, and are distinguished automatically:

    - `use_stellsym=True` (the usual case): `trgt_phi` covers *half* a field
      period and is offset by half a grid spacing. The data are first extended
      to a full field period using stellarator symmetry,
      `B_n(-theta, -phi) = -B_n(theta, phi)`.
    - `use_stellsym=False`: `trgt_phi` covers a full field period starting at
      phi=0, and is used as-is.
    """
    B = np.asarray(data["B_external_normal"], dtype=float)
    trgt_phi = data["trgt_phi"]
    trgt_theta = data["trgt_theta"]
    vc_nfp = data["nfp"]

    if vc_nfp != nfp:
        raise ValueError(
            f"The virtual-casing data has nfp={vc_nfp} but this surface has nfp={nfp}. "
            f"The virtual-casing calculation must be for the same plasma boundary."
        )
    if B.shape != (len(trgt_phi), len(trgt_theta)):
        raise ValueError(
            f"B_external_normal has shape {B.shape}, expected "
            f"{(len(trgt_phi), len(trgt_theta))} = (trgt_nphi, trgt_ntheta)."
        )

    theta_start = _uniform_grid_start(trgt_theta, 1.0, "theta")
    if abs(theta_start) > 1e-12:
        raise ValueError(
            f"The virtual-casing theta grid starts at {theta_start} rather than 0, "
            f"which is not a grid simsopt produces."
        )

    # Half a field period (stellarator-symmetric) or a whole one? For the half
    # period simsopt shifts the grid by half a spacing, so it never starts at 0.
    nphi = len(trgt_phi)
    # (the half-period grid starts at 0.25/(nfp*nphi), the full-period one at 0)
    half_period = trgt_phi[0] > 0.125 / (nfp * nphi)
    period = (0.5 if half_period else 1.0) / nfp
    phi_start = _uniform_grid_start(trgt_phi, period, "phi")

    if half_period:
        # B_n(theta_i, phi + 1/(2*nfp)) = -B_n(-theta_i, 1/(2*nfp) - phi). On the
        # index grid phi_j = (j + 1/2) / (2*nfp*nphi), the mirrored point of j is
        # nphi-1-j, and -theta_i is theta at index (-i) % ntheta.
        B = np.concatenate([B, -np.roll(B[::-1, ::-1], 1, axis=1)], axis=0)
        period *= 2

    # Fractional coordinates in which the (extended) samples sit at j/nphi, i/ntheta.
    u_targets = (zeta / (2 * np.pi) - phi_start) / period
    v_targets = theta / (2 * np.pi)
    # (nzeta, ntheta) -> REGCOIL's (ntheta, nzeta).
    return _spectral_interpolate(B, u_targets, v_targets).T


def bnormal_from_focus_modes(bfm, bfn, bfc, bfs, theta, zeta):
    """Evaluate B_normal from FOCUS-format Bn Fourier modes on a (theta, zeta) grid."""
    angle = bfm[:, None, None] * theta[None, :, None] - bfn[:, None, None] * zeta[None, None, :]
    return (bfc[:, None, None] * np.cos(angle) + bfs[:, None, None] * np.sin(angle)).sum(axis=0)
