"""`regcoil.cut`: cut discrete coils from a `Solution`'s total current
potential (Phase 10, ADR-029). Cutting is a computation that produces coil
geometry; plotting it (`regcoil.plot.plot_3d`/`coil_3d`) is downstream, and
the irreversible MAKEGRID file write stays out of the plot path. Replaces
`cutCoilsFromRegcoil` / `cut_saddle_coil`, plus the finite-thickness ribbon
geometry ported from `m20160811_01_plotCoilsFromRegcoil.m`.
"""

from __future__ import annotations

import numpy as np


class CutCoils:
    """`2 * coils_per_half_period * nfp` closed coil curves in Cartesian
    space, cut from one `Solution`'s total current potential.

    `curves` is a list of `(3, npoints)` arrays, each a closed curve (the
    first point repeated at the end). All coils carry the same `current`
    (`net_poloidal_current` split evenly across all coils). If `cut()` was
    called with `thickness`, `ribbons()` returns the finite-thickness corner
    geometry for each coil.
    """

    def __init__(self, curves, current, nfp, thickness=None):
        self.curves = curves
        self.current = current
        self.nfp = nfp
        self.thickness = thickness
        self._ribbons = None

    def __len__(self):
        return len(self.curves)

    def ribbons(self):
        """List of `(4, 3, npoints)` arrays: the four corner curves of the
        finite-thickness ribbon for each coil (offset by `thickness/2` along
        the local surface-normal and curve-binormal directions). Raises
        `ValueError` if `cut()` was called without `thickness`.
        """
        if self.thickness is None:
            raise ValueError("cut() was called without thickness; no ribbon geometry available")
        return self._ribbons

    def save_makegrid(self, path):
        """Write a MAKEGRID-format `coils.*` file (replaces
        `cutCoilsFromRegcoil`'s file output)."""
        with open(path, "w") as f:
            f.write(f"periods {self.nfp}\n")
            f.write("begin filament\n")
            f.write("mirror NIL\n")
            for curve in self.curves:
                n = curve.shape[1]
                for k in range(n - 1):
                    f.write(
                        f"{curve[0, k]:.14e} {curve[1, k]:.14e} {curve[2, k]:.14e} "
                        f"{self.current:.14e}\n"
                    )
                f.write(
                    f"{curve[0, -1]:.14e} {curve[1, -1]:.14e} {curve[2, -1]:.14e} "
                    "0.0 1 Modular\n"
                )
            f.write("end\n")

    def plot(self, fig=None):
        """Convenience delegate for `regcoil.plot.coil_3d(self, ...)`."""
        from . import plot

        return plot.coil_3d(self, fig=fig)


def _trace_contours(theta, zeta, Phi_norm, nfp, coils_per_half_period):
    """Trace `2*coils_per_half_period` closed contours of the (already
    period-normalized) potential `Phi_norm`, one per level, and replicate
    each across all `nfp` field periods by shifting `zeta`. Returns a list of
    `(theta_pts, zeta_pts)` pairs, one per coil (length
    `2*coils_per_half_period*nfp`).

    `Phi_norm` increases by exactly 1 across one field period (by
    construction, see `cut`); tracing on a 3-period-wide extension lets a
    contour that winds across the period seam close up correctly, matching
    the legacy `cutCoilsFromRegcoil` construction.
    """
    import contourpy

    d = 2 * np.pi / nfp
    zeta_3 = np.concatenate([zeta - d, zeta, zeta + d])
    Phi_3 = np.concatenate([Phi_norm - 1, Phi_norm, Phi_norm + 1], axis=1)

    ncoils_per_period = 2 * coils_per_half_period
    levels = np.linspace(0, 1, ncoils_per_period, endpoint=False)
    levels = levels + 0.5 * (levels[1] - levels[0])

    generator = contourpy.contour_generator(
        x=zeta_3, y=theta, z=Phi_3, line_type=contourpy.LineType.Separate,
    )

    curves_theta_zeta = []
    for level in levels:
        lines = generator.lines(float(level))
        if len(lines) != 1:
            raise RuntimeError(
                f"Expected exactly one contour at level {level!r} within one field "
                f"period, found {len(lines)}. The coil contours may cross the "
                "theta=0 line; try a different theta_shift."
            )
        curves_theta_zeta.append(lines[0])  # (npoints, 2): columns are (zeta, theta)

    surface_curves = []
    for zeta_theta in curves_theta_zeta:
        for jfp in range(nfp):
            zeta_pts = zeta_theta[:, 0] + d * jfp
            theta_pts = zeta_theta[:, 1]
            surface_curves.append((theta_pts, zeta_pts))
    return surface_curves


def _ribbon_corners(coil, theta_pts, zeta_pts, r_open, thickness):
    """Finite-thickness ribbon corners for one coil: offset the centerline
    `r_open` by `thickness/2` along the local surface-normal and
    curve-binormal directions (the Plotly port of the normal+binormal
    offset in `m20160811_01_plotCoilsFromRegcoil.m`)."""
    evaluated = coil.evaluate_at(theta_pts, zeta_pts)
    normal = np.cross(evaluated["drdzeta"], evaluated["drdtheta"], axis=0)
    normal /= np.linalg.norm(normal, axis=0)

    npts = r_open.shape[1]
    next_idx = np.roll(np.arange(npts), -1)
    prev_idx = np.roll(np.arange(npts), 1)
    tangent = r_open[:, next_idx] - r_open[:, prev_idx]
    tangent /= np.linalg.norm(tangent, axis=0)

    binormal = np.cross(tangent, normal, axis=0)
    binormal /= np.linalg.norm(binormal, axis=0)

    half = 0.5 * thickness
    corners = np.stack(
        [
            r_open + half * (normal + binormal),
            r_open + half * (normal - binormal),
            r_open + half * (-normal - binormal),
            r_open + half * (-normal + binormal),
        ]
    )  # (4, 3, npoints)
    return np.concatenate([corners, corners[:, :, :1]], axis=2)  # close the loop


def cut(solution, coils_per_half_period, thickness=None, theta_shift=0):
    """Cut `solution`'s total current potential into
    `2 * coils_per_half_period * nfp` discrete coils.

    Contours the total Phi at `2*coils_per_half_period` evenly-spaced levels
    (per field period), maps each contour's `(theta, zeta)` points to 3D via
    the coil surface's Fourier modes (`CoilSurface.evaluate_at`), and (if
    `thickness` is given) also builds finite-thickness ribbon geometry.
    Assumes the contours do not cross the theta=0 line after the
    `theta_shift` grid-point roll (matching the legacy `cutCoilsFromRegcoil`
    assumption -- adjust `theta_shift` if the default contours are rejected).

    Returns a `CutCoils`. `CutCoils.save_makegrid(path)` writes a MAKEGRID
    `coils.*` file.
    """
    prob = solution.problem
    coil = prob.coil
    nfp = coil.nfp
    net_poloidal_current = prob.net_poloidal_current

    theta = np.roll(coil.theta, theta_shift)
    theta = theta[0] + np.linspace(0, 2 * np.pi, len(theta), endpoint=False)
    zeta = coil.zeta

    Phi = np.roll(solution.current_potential(), theta_shift, axis=0)
    if abs(net_poloidal_current) > np.finfo(float).eps:
        Phi_norm = Phi / net_poloidal_current * nfp
    else:
        Phi_norm = Phi / np.max(np.abs(Phi))

    surface_curves = _trace_contours(theta, zeta, Phi_norm, nfp, coils_per_half_period)
    ncoils = len(surface_curves)
    current = net_poloidal_current / ncoils

    curves = []
    ribbons = [] if thickness is not None else None
    for theta_pts, zeta_pts in surface_curves:
        evaluated = coil.evaluate_at(theta_pts, zeta_pts)
        r_open = evaluated["r"]
        if ribbons is not None:
            ribbons.append(_ribbon_corners(coil, theta_pts, zeta_pts, r_open, thickness))
        curves.append(np.concatenate([r_open, r_open[:, :1]], axis=1))

    cut_coils = CutCoils(curves, current, nfp, thickness=thickness)
    cut_coils._ribbons = ribbons
    return cut_coils
