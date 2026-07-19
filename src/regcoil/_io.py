"""File readers for REGCOIL geometry/physics input formats.

Pure Python (no Fortran). VMEC `wout` reading tries
`scipy.io.netcdf_file` first (already a required dependency, and sufficient
for classic-format wout files) and falls back to the optional `netCDF4`
package for HDF5-based files.
"""

from __future__ import annotations

import numpy as np


def read_vmec_wout(filename):
    """Read the subset of a VMEC `wout` file that REGCOIL needs."""
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

    nfp = int(get("nfp"))
    if "lasym__logical__" in variables:
        lasym = bool(get("lasym__logical__"))
    else:
        lasym = bool(get("lasym"))

    xm = get("xm").astype(np.int64)
    xn = get("xn").astype(np.int64)
    rmnc = get("rmnc")
    zmns = get("zmns")
    rmns = get("rmns") if lasym and "rmns" in variables else None
    zmnc = get("zmnc") if lasym and "zmnc" in variables else None
    bvco = get("bvco")
    bsubvmnc = get("bsubvmnc")
    lmns = get("lmns") if "lmns" in variables else None
    mpol = int(get("mpol"))
    ntor = int(get("ntor"))

    return dict(
        nfp=nfp, lasym=lasym, xm=xm, xn=xn, rmnc=rmnc, zmns=zmns, rmns=rmns, zmnc=zmnc,
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


def bnormal_from_focus_modes(bfm, bfn, bfc, bfs, theta, zeta):
    """Evaluate B_normal from FOCUS-format Bn Fourier modes on a (theta, zeta) grid."""
    angle = bfm[:, None, None] * theta[None, :, None] - bfn[:, None, None] * zeta[None, None, :]
    return (bfc[:, None, None] * np.cos(angle) + bfs[:, None, None] * np.sin(angle)).sum(axis=0)
