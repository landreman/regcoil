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

# Typical usage

This page walks through a complete, realistic workflow: build surfaces, solve
a λ family, hit a target coil-current density, and save/load the result. Every
cell on this page is executed as part of the documentation build (ADR-030).

```{code-cell} ipython3
import numpy as np
import regcoil

ds = regcoil.examples("NCSX")
```

## 1. Surfaces

The plasma surface comes from a VMEC `wout` file; its target normal field
comes from a BNORM file (the contribution of in-plasma current to
`Bnormal`, which the coils must cancel):

```{code-cell} ipython3
plasma = regcoil.PlasmaSurface.from_wout(ds.wout, ntheta=24, nzeta=24)
plasma.set_bnormal_from_bnorm_file(ds.bnorm)
```

The coil ("winding") surface can be loaded independently from a NESCIN file
(`CoilSurface.from_nescin`), or computed as a uniform offset of the plasma
surface. The offset construction calls the Fortran kernel
(`regcoil_uniform_offset_surface`) once and returns ordinary Fourier
coefficients, so the result is a plain `CoilSurface` from then on:

```{code-cell} ipython3
coil = regcoil.CoilSurface.from_uniform_offset(
    plasma, separation=0.3, ntheta=24, nzeta=24, mpol=8, ntor=8
)
```

## 2. Assembling the problem

`Regcoil` builds the basis functions, calls the Fortran kernel
`build_g_and_h` once, and forms the eigendecomposition that makes every
subsequent λ solve or scan cheap:

```{code-cell} ipython3
problem = regcoil.Regcoil(plasma, coil, mpol_potential=8, ntor_potential=8)
```

## 3. Scanning λ

```{code-cell} ipython3
lambdas = np.array([0.0, 1e-16, 1e-15, 1e-14, 1e-13, 1e-12])
scan = problem.scan(lambdas)

for lam, f_B, max_K in zip(scan.lam, scan.f_B, scan.max_K):
    print(f"lam={lam:>8.1e}  f_B={f_B:.3e}  max|K|={max_K:.3e} A/m")
```

## 4. Hitting a target

Rather than scanning by hand, `solve_for_target` bisects in `log(lambda)` for
the λ whose `Solution.<metric>` equals a requested value — for example, the
λ that limits the peak coil-surface current density to a manufacturable
value:

```{code-cell} ipython3
target = problem.solve_for_target("max_K", 6.0e6)
print(f"lam = {target.lam:.3e}, max|K| = {target.max_K:.3e} A/m")
```

`solve_for_target` raises `ValueError` if the requested value is not between
the `lam=0` and `lam=inf` limits (i.e. the target is not achievable by any
choice of regularization).

## 5. Saving and loading

`regcoil.save`/`regcoil.load` write/read a self-contained NetCDF-4 file
(readable by generic tools like `xarray` or `h5py`, no `regcoil` import
required). Passing `solutions=scan` also carries along the `problem`,
`plasma`, and `coil` it depends on:

```{code-cell} ipython3
import tempfile, os

path = os.path.join(tempfile.mkdtemp(), "run.nc")
regcoil.save(path, plasma=plasma, coil=coil, problem=problem, solutions=scan)

result = regcoil.load(path)
print(len(result.solutions), "solutions loaded")
```

Loading is cheap and needs **no Fortran kernel and no expensive BLAS**: every
`Solution` already carries its own `current_potential`/`current_density`/
`Bnormal_total` grids, computed once at solve time and stored on disk.

```{code-cell} ipython3
sol = result.solutions[0]
print("current_potential shape:", sol.current_potential().shape)
print("current_density shape:", sol.current_density().shape)
```

A loaded `Regcoil` only rebuilds its operators (kernel call + eigendecomposition)
if you ask it to solve a *new* λ that wasn't in the saved scan.

## Next steps

See [Plotting](plotting.md) for every plot type `regcoil.plot` provides
(cross-sections, Pareto fronts, current potential/density maps, 3D views), and
[Cutting into discrete coils](plotting.md#cutting-into-discrete-coils) for
`regcoil.cut`.
