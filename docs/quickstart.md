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

# Quickstart

## The 30-second version

```{doctest}
>>> import regcoil
>>> ds = regcoil.examples("W7-X")
>>> plasma = regcoil.PlasmaSurface.from_wout(ds.wout, ntheta=64, nzeta=64)
>>> plasma.set_bnormal_from_bnorm_file(ds.bnorm)
>>> coil = regcoil.CoilSurface.from_uniform_offset(
...     plasma, separation=0.3, ntheta=64, nzeta=64, mpol=12, ntor=12
... )
>>> problem = regcoil.Regcoil(plasma, coil, mpol_potential=12, ntor_potential=12)
>>> solution = problem.solve(lam=1e-14)
>>> print(f"f_B = {solution.f_B:.1e}, f_K = {solution.f_K:.1e}")
f_B = 2.7e-01, f_K = 1.1e+15
```

The rest of this page walks through the same steps with more explanation.

## Bundled example data

`regcoil.examples` gives you ready-to-use VMEC (`wout`), plasma-current
(`bnorm`), and winding-surface (`nescin`) files, so you don't need external
data to try the package:

```{code-cell} ipython3
import regcoil

regcoil.examples.available()
```

```{code-cell} ipython3
ds = regcoil.examples("W7-X")
ds
```

## Building the plasma and coil surfaces

A {class}`~regcoil.PlasmaSurface` is built from a VMEC `wout` file, and its
target normal field from a BNORM file. A {class}`~regcoil.CoilSurface` can be
an independent NESCIN winding surface, or (as here) a uniform offset of the
plasma surface computed by the Fortran kernel:

```{code-cell} ipython3
plasma = regcoil.PlasmaSurface.from_wout(ds.wout, ntheta=64, nzeta=64)
plasma.set_bnormal_from_bnorm_file(ds.bnorm)

coil = regcoil.CoilSurface.from_uniform_offset(
    plasma, separation=0.3, ntheta=64, nzeta=64, mpol=12, ntor=12
)
print(f"plasma nfp={plasma.nfp}, plasma area={plasma.area:.3g} m^2")
print(f"coil area={coil.area:.3g} m^2")
```

(These resolutions are kept low so this page builds quickly; see
[Usage](usage.md) and [Plotting](plotting.md) for more realistic settings.)

## Solving

A {class}`~regcoil.Regcoil` problem assembles the inductance matrices and
their eigendecomposition once; solving for a particular regularization
parameter λ, or scanning over many λ values, is then cheap:

```{code-cell} ipython3
problem = regcoil.Regcoil(plasma, coil, mpol_potential=12, ntor_potential=12)
solution = problem.solve(lam=1e-14)
print(f"f_B = {solution.f_B:.3e}")
print(f"f_K = {solution.f_K:.3e}")
print(f"max|K| = {solution.max_K:.3e} A/m")
print(f"rms|K| = {solution.rms_K:.3e} A/m")
print(f"max|Bnormal| = {solution.max_Bnormal:.3e} T")
print(f"max|Bnormal|/|B| = {solution.max_Bnormal_over_B:.3e}")
print(f"avg|Bnormal|/|B| = {solution.avg_Bnormal_over_B:.3e}")
```

```{code-cell} ipython3
import numpy as np

lambdas = np.array([0.0, 1e-16, 1e-15, 1e-14, 1e-12])
scan = problem.scan(lambdas)
for lam, f_B, f_K in zip(scan.lam, scan.f_B, scan.f_K):
    print(f"lam={lam:.1e}  f_B={f_B:.3e}  f_K={f_K:.3e}")
```

Smaller λ favors a smaller field residual (`f_B`) at the cost of a more
complex current potential (larger `f_K`); this trade-off is the Pareto front
plotted in [Plotting](plotting.md#pareto-fronts).

## Next steps

- [Usage](usage.md): a fuller workflow, including `solve_for_target` and
  saving/loading a run.
- [Plotting](plotting.md): every plot type `regcoil.plot` provides.
- [API reference](api.md): the full object model.
