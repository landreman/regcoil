---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.19.4
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---

# Plotting

`regcoil.plot` has one atomic function per plot type; every function accepts
`ax=` (matplotlib) or `fig=` (Plotly) and returns it, so panels compose into
larger dashboards with no duplicated logic. Each `obj.plot_*()` convenience
method (e.g. `Solution.plot_current_potential()`) delegates to the matching
free function. matplotlib/Plotly are imported lazily inside these functions,
so `import regcoil` never pulls in a plotting library.

```{code-cell} ipython3
%matplotlib inline

import numpy as np
import plotly.io as pio
import regcoil

# Self-contained <script>/<div> HTML output (loads plotly.js from a CDN in
# the reader's browser) instead of the ipywidget-style mimebundle Jupyter's
# own frontend understands but a static Sphinx page cannot render.
pio.renderers.default = "notebook_connected"

ds = regcoil.examples("W7-X")
plasma = regcoil.PlasmaSurface.from_wout(ds.wout, ntheta=64, nzeta=64)
plasma.set_bnormal_from_bnorm_file(ds.bnorm)
coil = regcoil.CoilSurface.from_uniform_offset(
    plasma, separation=0.3, ntheta=64, nzeta=64, mpol=12, ntor=12
)
problem = regcoil.Regcoil(plasma, coil, mpol_potential=12, ntor_potential=12)

lambdas = np.array([0.0, 1e-16, 1e-15, 1e-14, 1e-13, 1e-12])
scan = problem.scan(lambdas)
solution = scan[3]  # one representative lambda for the single-solution plots
```

## Surface cross sections

`cross_sections` shows the plasma and coil surface together, at several
toroidal angles; `cross_sections_overlay` shows one surface's cross sections
overlaid and colored by angle. Both go through `Surface.cross_section(phi)`,
which is correct regardless of the surface's `standard_toroidal_angle`
convention (it derives `phi = atan2(y, x)` from the real Cartesian grid
instead of assuming constant-ζ means constant physical angle).

```{code-cell} ipython3
regcoil.plot.cross_sections(plasma, coil);
```

```{code-cell} ipython3
regcoil.plot.cross_sections_overlay(plasma);
```

## Pareto fronts

Trade-off between field error (`f_B`) and current-potential complexity
(`f_K`) across a λ scan. `pareto` accepts a single `SolutionScan` or a list of
them (overlaid, e.g. to compare two coil separations):

```{code-cell} ipython3
regcoil.plot.pareto([scan]);
```

Any pair of `f_B`, `f_K`, `max_K`, `max_Bnormal` can be plotted against each
other via `x=`/`y=`.

## f_B / f_K vs. lambda

```{code-cell} ipython3
regcoil.plot.lambda_scan(scan);
```

## Current potential

`(theta, zeta)` contour map of the current potential Phi, for one `Solution`.
`kind='single_valued'` (default) excludes the secular net-current term;
`kind='total'` includes it:

```{code-cell} ipython3
regcoil.plot.current_potential(solution, kind="total");
```

`current_potential_scan` lays out a panel grid across (up to `nmax`) lambdas
in one scan:

```{code-cell} ipython3
regcoil.plot.current_potential_scan(scan, kind="total");
```

## Current density

```{code-cell} ipython3
regcoil.plot.current_density(solution);
```

## B_normal

`component` is `'plasma_current'`, `'net_coil'`, or `'total'` (default).

```{code-cell} ipython3
regcoil.plot.bnormal(solution, component="total");
```

```{code-cell} ipython3
regcoil.plot.bnormal_scan(scan);
```

## 3D interactive views

`plot_3d` renders any subset of {plasma surface, winding surface, cut coils}
in one Plotly scene; `winding_surface_style` is `'wireframe'` (default),
`'translucent'`, or `'solid'`:

```{code-cell} ipython3
regcoil.plot.plot_3d(plasma=plasma, winding_surface=coil, winding_surface_style="translucent")
```

## Cutting into discrete coils

`regcoil.cut.cut` contours the total current potential into discrete coils
and (optionally) finite-thickness ribbon geometry; cutting is a separate,
explicit step from plotting, and its `CutCoils.save_makegrid` is the only
irreversible file write on this page:

```{code-cell} ipython3
cut_coils = regcoil.cut.cut(solution, coils_per_half_period=5, thickness=0.05)
regcoil.plot.coil_3d(cut_coils)
```

```{code-cell} ipython3
regcoil.plot.plot_3d(plasma=plasma, coils=cut_coils)
```

```{code-cell} ipython3
import tempfile, os

path = os.path.join(tempfile.mkdtemp(), "coils.example")
cut_coils.save_makegrid(path)
print(f"wrote {os.path.getsize(path)} bytes to a MAKEGRID coils file")
```

## Everything at once

`regcoil.plot.all` (alias `dashboard`) composes every atomic function above
from anything exposing `.plasma`, `.coil`, `.solutions` (any may be absent) —
a live `Regcoil` doesn't carry a `.solutions` attribute itself, but the
container returned by `regcoil.load()` does, so it's the natural input here:

```{code-cell} ipython3
import tempfile, os

path = os.path.join(tempfile.mkdtemp(), "run.nc")
regcoil.save(path, plasma=plasma, coil=coil, problem=problem, solutions=scan)
result = regcoil.load(path)

figures = regcoil.plot.all(result)
print(f"{len(figures)} figures")
```

## The `regcoil` CLI

The same plots are reachable without writing any Python, via the console
script installed with the package (`regcoil plot|compare|cut`):

```bash
regcoil plot run.nc                      # full dashboard
regcoil plot run.nc --only pareto        # one atomic plot
regcoil compare run_a.nc run_b.nc        # overlay Pareto fronts across runs
regcoil cut run.nc --coils-per-half-period 3 --thickness 0.05 --out coils.out
```

Run `regcoil --help` (or `regcoil plot --help`, etc.) for the full option
list.
