"""`regcoil` command-line interface: one console script
with subcommands, each a thin `load() -> plot/cut -> show/save` wrapper. File
paths are the CLI's concern only -- the library functions stay object-only.
Replaces the separate `regcoilPlot`/`compareRegcoil`/`cutCoilsFromRegcoil`/
`cut_saddle_coil` scripts.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

ATOMIC_PLOTS = (
    "cross_sections", "cross_sections_overlay", "pareto", "lambda_scan",
    "current_potential", "current_density", "bnormal", "plot_3d",
)


def _load(path):
    from . import _serialize
    return _serialize.load(path)


def _is_plotly_figure(obj):
    return type(obj).__module__.startswith("plotly")


def _dispatch_figures(figs, save, no_show):
    import matplotlib.pyplot as plt

    mpl_figs = []
    plotly_figs = []
    for f in figs:
        if _is_plotly_figure(f):
            plotly_figs.append(f)
        elif hasattr(f, "savefig"):
            mpl_figs.append(f)  # already a Figure
        else:
            mpl_figs.append(f.figure)  # an Axes -> its parent Figure

    if save:
        total = len(mpl_figs) + len(plotly_figs)
        if total == 1:
            if plotly_figs:
                plotly_figs[0].write_html(save)
            else:
                mpl_figs[0].savefig(save)
            print(f"Wrote {save}")
        else:
            root, dot, ext = save.rpartition(".")
            root = root if dot else save
            for i, fig in enumerate(mpl_figs):
                out = f"{root}_{i}.pdf"
                fig.savefig(out)
                print(f"Wrote {out}")
            for i, fig in enumerate(plotly_figs):
                out = f"{root}_plotly{i}.html"
                fig.write_html(out)
                print(f"Wrote {out}")

    if not no_show:
        for fig in plotly_figs:
            fig.show()
        if mpl_figs:
            plt.show()


def _cmd_plot(args):
    from . import plot

    data = _load(args.file)

    if not args.only:
        figs = plot.all(data)
        _dispatch_figures(figs, args.save, args.no_show)
        return

    figs = []
    for name in args.only:
        if name == "cross_sections":
            phi = None if args.phi is None else np.asarray(args.phi, dtype=float)
            figs.append(plot.cross_sections(data.plasma, data.coil, phi=phi))
        elif name == "cross_sections_overlay":
            phi = None if args.phi is None else np.asarray(args.phi, dtype=float)
            surf = data.plasma if args.surface == "plasma" else data.coil
            figs.append(plot.cross_sections_overlay(surf, phi=phi).figure)
        elif name == "pareto":
            figs.append(plot.pareto(data.solutions).figure)
        elif name == "lambda_scan":
            figs.append(plot.lambda_scan(data.solutions).figure)
        elif name in ("current_potential", "current_density", "bnormal"):
            sol = data.solutions[args.lambda_index]
            fn = getattr(plot, name)
            figs.append(fn(sol).figure)
        elif name == "plot_3d":
            figs.append(plot.plot_3d(plasma=data.plasma, winding_surface=data.coil))
        else:  # pragma: no cover - argparse `choices` already restricts this
            raise SystemExit(f"Unknown --only {name!r}; choose from {ATOMIC_PLOTS}")

    _dispatch_figures(figs, args.save, args.no_show)


def _cmd_compare(args):
    from . import plot

    scans = [_load(f).solutions for f in args.files]
    ax = plot.pareto(scans, x=args.x, y=args.y, labels=args.files)
    _dispatch_figures([ax.figure], args.save, args.no_show)


def _cmd_cut(args):
    from . import cut as cut_module
    from . import plot

    data = _load(args.file)
    sol = data.solutions[args.lambda_index]
    cut_coils = cut_module.cut(
        sol, args.coils_per_half_period, thickness=args.thickness, theta_shift=args.theta_shift
    )

    out = args.out or f"coils.{Path(args.file).stem}"
    cut_coils.save_makegrid(out)
    print(f"Wrote {out} ({len(cut_coils)} coils, {cut_coils.current:.6g} A each)")

    if args.plot:
        fig = plot.plot_3d(winding_surface=data.coil, coils=cut_coils)
        _dispatch_figures([fig], None, args.no_show)


def build_parser():
    parser = argparse.ArgumentParser(prog="regcoil")
    sub = parser.add_subparsers(dest="command", required=True)

    p_plot = sub.add_parser("plot", help="Plot a saved run (regcoil.plot.all, or --only atomic plots)")
    p_plot.add_argument("file")
    p_plot.add_argument(
        "--only", action="append", choices=ATOMIC_PLOTS,
        help="Plot only this atomic function (repeatable); default is the full dashboard",
    )
    p_plot.add_argument("--phi", type=float, nargs="+", help="Cross-section physical angles (radians)")
    p_plot.add_argument(
        "--surface", default="plasma", choices=("plasma", "coil"),
        help="Which surface to plot for --only cross_sections_overlay (default: plasma)",
    )
    p_plot.add_argument("--lambda-index", type=int, default=0, help="Solution index for the (theta,zeta) maps")
    p_plot.add_argument("--save", help="Save to OUT.pdf (matplotlib) or OUT.html (Plotly)")
    p_plot.add_argument("--no-show", action="store_true", help="Skip the interactive window")
    p_plot.set_defaults(func=_cmd_plot)

    p_compare = sub.add_parser("compare", help="Overlay Pareto fronts across several saved runs")
    p_compare.add_argument("files", nargs="+")
    p_compare.add_argument("--x", default="f_K", choices=("f_B", "f_K", "max_K", "max_Bnormal"))
    p_compare.add_argument("--y", default="f_B", choices=("f_B", "f_K", "max_K", "max_Bnormal"))
    p_compare.add_argument("--save", help="Save to OUT.pdf (matplotlib) or OUT.html (Plotly)")
    p_compare.add_argument("--no-show", action="store_true")
    p_compare.set_defaults(func=_cmd_compare)

    p_cut = sub.add_parser("cut", help="Cut coils from a saved run and write a MAKEGRID coils.* file")
    p_cut.add_argument("file")
    p_cut.add_argument("--coils-per-half-period", type=int, required=True)
    p_cut.add_argument("--lambda-index", type=int, default=0)
    p_cut.add_argument("--thickness", type=float, default=None)
    p_cut.add_argument("--theta-shift", type=int, default=0)
    p_cut.add_argument("--out", help="Output coils.* filename (default: coils.<FILE stem>)")
    p_cut.add_argument("--save", help="Save to OUT.pdf (matplotlib) or OUT.html (Plotly)")
    p_cut.add_argument("--plot", action="store_true", help="Also show the 3D finite-thickness coil figure")
    p_cut.add_argument("--no-show", action="store_true")
    p_cut.set_defaults(func=_cmd_cut)

    return parser


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    sys.exit(main())
