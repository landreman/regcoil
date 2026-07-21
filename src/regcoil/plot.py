"""`regcoil.plot`: visualization for the object model (Phase 10, ADR-029).

Consumes only in-memory objects (`PlasmaSurface`/`CoilSurface`/`Regcoil`/
`Solution`/`SolutionScan`, or the `regcoil.load()` container) -- never file
paths. matplotlib/plotly are imported lazily inside each function, so
`import regcoil` stays plotting-free. Every atomic function accepts `ax=`
(matplotlib) or `fig=` (plotly) and returns it; the multi-panel/scan
functions and `all()` are pure composition of the atomic functions (no
plotting logic of their own).
"""

from __future__ import annotations

import numpy as np

from .regcoil import _unflatten_grid

#: Default matplotlib figure size, used whenever a function must create its
#: own figure (i.e. no `ax=` was supplied).
DEFAULT_FIGSIZE = (14.5, 8.1)


def _new_ax(ax):
    import matplotlib.pyplot as plt

    if ax is None:
        _, ax = plt.subplots(figsize=DEFAULT_FIGSIZE)
    return ax


def _make_grid_axes(fig, nrows, ncols):
    return np.array(
        [fig.add_subplot(nrows, ncols, i + 1) for i in range(nrows * ncols)]
    ).reshape(nrows, ncols)


def _grid_shape(n):
    ncols = int(np.ceil(np.sqrt(n)))
    nrows = int(np.ceil(n / ncols))
    return nrows, ncols


def _pick_lambda_indices(solutions, nmax):
    """Indices into `solutions`, sorted by lambda and subsampled to at most
    `nmax` entries evenly spaced along the sorted sequence (matching the
    legacy `regcoilPlot`'s ilambda_to_plot selection)."""
    order = np.argsort(solutions.lam)
    n = len(order)
    if n <= nmax:
        return order
    pick = np.unique(np.round(np.linspace(0, n - 1, nmax)).astype(int))
    return order[pick]


def cross_sections_overlay(surf, phi=None, ax=None):
    """One surface's cross section(s) at fixed physical toroidal angle(s),
    overlaid on a single axes and colored by `phi`, via `Surface.cross_section`
    (correct for both `standard_toroidal_angle` values, ADR-025). Default
    `phi = [0, 0.5, 1, 1.5] * pi / nfp`.
    """
    ax = _new_ax(ax)
    if phi is None:
        phi = np.array([0, 0.5, 1, 1.5]) * np.pi / surf.nfp
    phi = np.atleast_1d(np.asarray(phi, dtype=float))

    import matplotlib.pyplot as plt

    colors = plt.cm.viridis(np.linspace(0, 0.85, len(phi)))

    R, Z = surf.cross_section(phi)
    for i in range(len(phi)):
        r_closed = np.append(R[i], R[i, 0])
        z_closed = np.append(Z[i], Z[i, 0])
        ax.plot(r_closed, z_closed, ".-", color=colors[i], label=f"phi={phi[i]:.3f}")

    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("R [m]")
    ax.set_ylabel("Z [m]")
    ax.legend(fontsize="x-small")
    return ax


def cross_sections(plasma, coil, phi=None, fig=None):
    """Grid of subplots, one per physical toroidal angle in `phi`, each
    showing the plasma (red) and coil (blue) cross section at that angle
    (composition of `cross_sections_overlay`'s underlying geometry, but with
    surface-identity color coding instead of phi color coding, to avoid the
    clutter of overlaying every phi on one axes). Default
    `phi = [0, 0.5, 1, 1.5] * pi / nfp`.
    """
    import matplotlib.pyplot as plt

    if phi is None:
        phi = np.array([0, 0.5, 1, 1.5]) * np.pi / plasma.nfp
    phi = np.atleast_1d(np.asarray(phi, dtype=float))

    R_plasma, Z_plasma = plasma.cross_section(phi)
    R_coil, Z_coil = coil.cross_section(phi)

    nrows, ncols = _grid_shape(len(phi))
    if fig is None:
        fig = plt.figure(figsize=DEFAULT_FIGSIZE)
    axes = _make_grid_axes(fig, nrows, ncols)
    for ax, i in zip(axes.flat, range(len(phi))):
        r_p = np.append(R_plasma[i], R_plasma[i, 0])
        z_p = np.append(Z_plasma[i], Z_plasma[i, 0])
        r_c = np.append(R_coil[i], R_coil[i, 0])
        z_c = np.append(Z_coil[i], Z_coil[i, 0])
        ax.plot(r_p, z_p, ".-", color="red", label="plasma")
        ax.plot(r_c, z_c, ".-", color="blue", label="coil")
        ax.set_aspect("equal", adjustable="box")
        ax.set_xlabel("R [m]")
        ax.set_ylabel("Z [m]")
        ax.set_title(f"phi={phi[i]:.3f}", fontsize="small")
        ax.legend(fontsize="x-small")
    fig.tight_layout()
    return fig


def pareto(scans, x="f_K", y="f_B", labels=None, ax=None):
    """Pareto front(s): `x`/`y` chosen among `f_B`, `f_K`, `max_K`,
    `max_Bnormal`. `scans` is a single `SolutionScan` or a list of them
    (overlaid, one line per scan)."""
    ax = _new_ax(ax)
    if not isinstance(scans, (list, tuple)):
        scans = [scans]
    for i, scan in enumerate(scans):
        label = labels[i] if labels is not None else None
        xv = np.asarray(getattr(scan, x))
        yv = np.asarray(getattr(scan, y))
        order = np.argsort(xv)
        ax.loglog(xv[order], yv[order], ".-", label=label)
    ax.set_xlabel(x)
    ax.set_ylabel(y)
    ax.grid(True, which="both", alpha=0.3)
    if labels is not None:
        ax.legend(fontsize="small")
    return ax


def lambda_scan(solutions, ax=None):
    """`f_B`/`f_K` vs lambda traces (the rest of the legacy scan grid), on
    twin y-axes of one shared `ax`."""
    ax = _new_ax(ax)
    lam = np.asarray(solutions.lam)
    order = np.argsort(lam)
    lam_sorted = lam[order]

    line_B, = ax.loglog(lam_sorted, np.asarray(solutions.f_B)[order], "C0.-", label="f_B")
    ax.set_xlabel("lambda")
    ax.set_ylabel("f_B", color="C0")
    ax.tick_params(axis="y", colors="C0")

    ax2 = ax.twinx()
    line_K, = ax2.loglog(lam_sorted, np.asarray(solutions.f_K)[order], "C1.-", label="f_K")
    ax2.set_ylabel("f_K", color="C1")
    ax2.tick_params(axis="y", colors="C1")

    ax.legend(handles=[line_B, line_K], fontsize="small")
    return ax


def _single_valued_current_potential(solution):
    prob = solution.problem
    coil = prob.coil
    return _unflatten_grid(prob.basis_functions @ solution.solution, coil.ntheta, coil.nzeta)


def current_potential(solution, kind="single_valued", ax=None):
    """`(theta, zeta)` contour map of the current potential Phi, one
    `Solution` (one lambda) at a time. `kind` is `'single_valued'` (default)
    or `'total'` (includes the secular net-current term)."""
    ax = _new_ax(ax)
    coil = solution.problem.coil
    if kind == "single_valued":
        data = _single_valued_current_potential(solution)
    elif kind == "total":
        data = solution.current_potential()
    else:
        raise ValueError(f"kind must be 'single_valued' or 'total', got {kind!r}")

    cf = ax.contourf(coil.zeta, coil.theta, data, 20)
    ax.figure.colorbar(cf, ax=ax)
    ax.set_xlabel("zeta")
    ax.set_ylabel("theta")
    ax.set_title(f"current potential ({kind}), lambda={solution.lam:.3g}", fontsize="small")
    return ax


def current_density(solution, ax=None):
    """`(theta, zeta)` contour map of `|K|`, one `Solution` at a time."""
    ax = _new_ax(ax)
    coil = solution.problem.coil
    K = solution.current_density()
    magnitude = np.sqrt(np.sum(K * K, axis=0))
    cf = ax.contourf(coil.zeta, coil.theta, magnitude, 20)
    ax.figure.colorbar(cf, ax=ax)
    ax.set_xlabel("zeta")
    ax.set_ylabel("theta")
    ax.set_title(f"|K| [A/m], lambda={solution.lam:.3g}", fontsize="small")
    return ax


def bnormal(solution, component="total", ax=None):
    """`(theta, zeta)` contour map of B_normal on the plasma surface.
    `component` is `'plasma_current'`, `'net_coil'`, or `'total'` (default).
    `'net_coil'` needs no Fortran kernel on a loaded run (ADR-028/ADR-029:
    `Regcoil.Bnormal_from_net_coil_currents` is stored on disk).
    """
    ax = _new_ax(ax)
    prob = solution.problem
    plasma = prob.plasma
    if component == "plasma_current":
        data = plasma.Bnormal_from_plasma_current
    elif component == "net_coil":
        data = prob.Bnormal_from_net_coil_currents
    elif component == "total":
        data = solution.Bnormal_total
    else:
        raise ValueError(
            f"component must be 'plasma_current', 'net_coil', or 'total', got {component!r}"
        )

    cf = ax.contourf(plasma.zeta, plasma.theta, data, 20)
    ax.figure.colorbar(cf, ax=ax)
    ax.set_xlabel("zeta")
    ax.set_ylabel("theta")
    title = f"Bnormal ({component})"
    if component == "total":
        title += f", lambda={solution.lam:.3g}"
    ax.set_title(title, fontsize="small")
    return ax


def current_potential_scan(solutions, kind="total", nmax=16, fig=None):
    """Multi-lambda panel grid of `current_potential` (composition of the
    atomic function): up to `nmax` lambdas, evenly spaced by rank."""
    import matplotlib.pyplot as plt

    idx = _pick_lambda_indices(solutions, nmax)
    nrows, ncols = _grid_shape(len(idx))
    if fig is None:
        fig = plt.figure(figsize=DEFAULT_FIGSIZE)
    axes = _make_grid_axes(fig, nrows, ncols)
    for ax, i in zip(axes.flat, idx):
        current_potential(solutions[i], kind=kind, ax=ax)
    fig.suptitle(f"current potential ({kind})")
    fig.tight_layout()
    return fig


def bnormal_scan(solutions, nmax=16, fig=None):
    """Multi-lambda panel grid of `bnormal` (composition of the atomic
    function): the plasma-current and net-coil-current panels, plus up to
    `nmax` total-Bnormal panels evenly spaced by lambda rank."""
    import matplotlib.pyplot as plt

    idx = _pick_lambda_indices(solutions, nmax)
    n_panels = len(idx) + 2
    nrows, ncols = _grid_shape(n_panels)
    if fig is None:
        fig = plt.figure(figsize=DEFAULT_FIGSIZE)
    axes = _make_grid_axes(fig, nrows, ncols)
    flat = axes.flat
    bnormal(solutions[0], component="plasma_current", ax=next(flat))
    bnormal(solutions[0], component="net_coil", ax=next(flat))
    for ax, i in zip(flat, idx):
        bnormal(solutions[i], component="total", ax=ax)
    fig.suptitle("Bnormal [T]")
    fig.tight_layout()
    return fig


def _close_surface_grid(surf):
    r = surf.r
    r = np.concatenate([r, r[:, :1, :]], axis=1)
    r = np.concatenate([r, r[:, :, :1]], axis=2)
    return r


def _add_surface_trace(fig, surf, style, color, name):
    import plotly.graph_objects as go

    r = _close_surface_grid(surf)
    X, Y, Z = r[0], r[1], r[2]

    if style == "wireframe":
        # Subsample rings so the wireframe stays legible.
        ntheta, nzetal = X.shape
        theta_step = max(1, ntheta // 30)
        zeta_step = max(1, nzetal // 60)
        for i in range(0, ntheta, theta_step):
            fig.add_trace(
                go.Scatter3d(
                    x=X[i, :], y=Y[i, :], z=Z[i, :], mode="lines",
                    line=dict(color=color, width=2), name=name, showlegend=False,
                )
            )
        for j in range(0, nzetal, zeta_step):
            fig.add_trace(
                go.Scatter3d(
                    x=X[:, j], y=Y[:, j], z=Z[:, j], mode="lines",
                    line=dict(color=color, width=2), name=name, showlegend=False,
                )
            )
    elif style in ("translucent", "solid"):
        opacity = 0.4 if style == "translucent" else 1.0
        fig.add_trace(
            go.Surface(
                x=X, y=Y, z=Z, showscale=False, opacity=opacity,
                colorscale=[[0, color], [1, color]], name=name,
            )
        )
    else:
        raise ValueError(
            f"style must be 'wireframe', 'translucent', or 'solid', got {style!r}"
        )


def _add_coils_trace(fig, coils):
    import plotly.graph_objects as go

    for curve in coils.curves:
        fig.add_trace(
            go.Scatter3d(
                x=curve[0], y=curve[1], z=curve[2], mode="lines",
                line=dict(color="black", width=4), showlegend=False,
            )
        )
    if coils.thickness is not None:
        for corners in coils.ribbons():
            X, Y, Z = corners[:, 0, :], corners[:, 1, :], corners[:, 2, :]
            fig.add_trace(
                go.Surface(
                    x=X, y=Y, z=Z, showscale=False, opacity=1.0,
                    colorscale=[[0, "orange"], [1, "orange"]],
                )
            )


def plot_3d(plasma=None, winding_surface=None, coils=None, winding_surface_style="wireframe", fig=None):
    """One interactive 3D Plotly scene, plotting any subset of {`plasma`
    surface, `winding_surface`, cut `coils`} -- omit what you don't want.
    `winding_surface_style` is `'wireframe'` (default), `'translucent'`, or
    `'solid'`.
    """
    import plotly.graph_objects as go

    if fig is None:
        fig = go.Figure()
    if plasma is not None:
        _add_surface_trace(fig, plasma, "solid", "hotpink", "plasma")
        _add_surface_trace(fig, plasma, "wireframe", "darkred", "plasma")
    if winding_surface is not None:
        _add_surface_trace(fig, winding_surface, winding_surface_style, "steelblue", "winding surface")
    if coils is not None:
        _add_coils_trace(fig, coils)
    fig.update_layout(scene=dict(aspectmode="data"))
    return fig


def coil_3d(coils, fig=None):
    """3D Plotly view of cut coils only (a convenience wrapper around
    `plot_3d(coils=coils)`)."""
    return plot_3d(coils=coils, fig=fig)


def all(data, nmax=16):
    """Everything: pure composition of the atomic functions (replaces the
    legacy `regcoilPlot` dashboard). `data` is anything exposing `.plasma`,
    `.coil`, `.solutions` (any may be `None`) -- e.g. a `regcoil.load()`
    container. Returns a list of figures: matplotlib figures for the 2D
    panels, plus one Plotly 3D figure if a plasma or coil surface is
    present.
    """
    import matplotlib.pyplot as plt

    figs = []
    plasma = getattr(data, "plasma", None)
    coil = getattr(data, "coil", None)
    solutions = getattr(data, "solutions", None)

    if plasma is not None:
        phis = np.array([0, 0.5, 1, 1.5]) * np.pi / plasma.nfp
        if coil is not None:
            fig = cross_sections(plasma, coil, phi=phis)
        else:
            fig = cross_sections_overlay(plasma, phi=phis).figure
        figs.append(fig)

    if solutions is not None:
        fig, ax = plt.subplots(figsize=DEFAULT_FIGSIZE)
        pareto(solutions, ax=ax)
        figs.append(fig)

        fig, ax = plt.subplots(figsize=DEFAULT_FIGSIZE)
        lambda_scan(solutions, ax=ax)
        figs.append(fig)

        figs.append(current_potential_scan(solutions, kind="total", nmax=nmax))
        figs.append(bnormal_scan(solutions, nmax=nmax))

    if plasma is not None or coil is not None:
        figs.append(plot_3d(plasma=plasma, winding_surface=coil))

    return figs


dashboard = all
