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
a λ family, hit a target coil-current density, and save/load the result.

```{code-cell} ipython3
import numpy as np
import regcoil

ds = regcoil.examples("W7-X")
```

## 1. Surfaces

The plasma surface comes from a VMEC `wout` file; its target normal field
comes from a virtual casing or bnorm file (the contribution of in-plasma current to
`Bnormal`, which the coils must cancel):

```{code-cell} ipython3
plasma = regcoil.PlasmaSurface.from_wout(ds.wout, ntheta=64, nzeta=64)
plasma.set_bnormal_from_virtual_casing(ds.vcasing)
# or plasma.set_bnormal_from_bnorm_file(ds.bnorm)
```

If you computed the plasma's normal field with
[simsopt](https://github.com/hiddenSymmetries/simsopt)'s virtual-casing module,
use `set_bnormal_from_virtual_casing`.
This command accepts either the path to a `vcasing*.nc`
file or a `simsopt.mhd.VirtualCasing` object, and works whether the
calculation was run with `use_stellsym=True` (the usual case) or `False`:

```python
plasma.set_bnormal_from_virtual_casing("vcasing_li383.nc")
```

Files are read with NetCDF directly, so simsopt need not be installed. Note
that simsopt gives `B_external_normal` on its own fixed grid rather than as
Fourier modes, so unless that grid happens to coincide with the plasma
surface's, the data are interpolated onto `ntheta` x `nzeta`. The
interpolation is trigonometric and so converges spectrally, but it cannot
invent structure the virtual-casing grid did not resolve: choose
`trgt_nphi`/`trgt_ntheta` in the virtual-casing run with the resolution you
want here in mind.

The coil ("winding") surface can be loaded independently from a NESCIN file
(`CoilSurface.from_nescin`), or computed as a uniform offset of the plasma
surface. The offset construction samples the offset points once and returns
ordinary Fourier coefficients, so the result is a plain `CoilSurface` from
then on:

```{code-cell} ipython3
coil = regcoil.CoilSurface.from_uniform_offset(
    plasma, separation=0.3, ntheta=64, nzeta=64, mpol=12, ntor=12
)
```

### Reparameterizing the poloidal angle

If the offset surface simply inherited its poloidal angle from the plasma,
that angle would often be poorly distributed -- points bunch up where the
offset compresses the surface and spread out elsewhere. `from_uniform_offset`
therefore relabels the points by default (`theta_reparameterization=
"uniform_arclength"`), leaving the shape untouched but making the angle better
behaved. Pass `theta_reparameterization=None` to keep the plasma's angle:

```{code-cell} ipython3
import numpy as np

inherited = regcoil.CoilSurface.from_uniform_offset(
    plasma, separation=0.3, ntheta=64, nzeta=64, mpol=12, ntor=12,
    theta_reparameterization=None,
)

def arclength_spread(surf):
    dr = surf.drdtheta[:, :, : surf.nzeta]
    speed = np.sqrt(np.sum(dr * dr, axis=0))
    return ((speed.max(axis=0) - speed.min(axis=0)) / speed.mean(axis=0)).max()

print(f"inherited angle:         {arclength_spread(inherited):.3f}")
print(f"uniform-arclength angle: {arclength_spread(coil):.3f}")
```

`"uniform_arclength"` makes `|dr/dtheta|` independent of `theta`;
`"curvature"` instead concentrates points where the cross section bends
sharply. For finer control pass a
{class}`~regcoil.UniformArclength` or {class}`~regcoil.CurvatureWeighted`
instance. The resulting {class}`~regcoil.ThetaMap` is kept on the surface as
`coil.theta_map`, and its `diagnostics` report how well the target was met.

Because the current-potential basis `sin(m*theta - n*zeta)` is defined in terms of the
coil surface's `theta`, this reparameterization changes the solution space.

## 2. Assembling the problem

`Regcoil` builds the basis functions, calls the Fortran kernel
`build_g_and_h` once, and forms the eigendecomposition that makes every
subsequent λ solve or scan cheap:

```{code-cell} ipython3
problem = regcoil.Regcoil(plasma, coil, mpol_potential=12, ntor_potential=12)
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
target = problem.solve_for_target("max_K", 4.0e6)
print(f"lam = {target.lam:.3e}, max|K| = {target.max_K:.3e} A/m")
```

Any of these `Solution` scalars can be used as the target metric:

- `f_B` — ∫ (B·n)² over the plasma surface
- `f_K` — ∫ K² over the coil surface
- `max_K` — peak |K| on the coil surface
- `rms_K` — RMS |K| on the coil surface
- `max_Bnormal` — max |B·n| on the plasma surface
- `max_Bnormal_over_B` — max( |B·n| / |B| ) on the plasma surface
- `avg_Bnormal_over_B` — area-mean( |B·n| / |B| ) on the plasma surface

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

Loading is cheap: every
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
