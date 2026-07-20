# Target API sketch

Design target for the migration—not a frozen public API. Refine via [DECISIONS.md](DECISIONS.md).

> **Architecture pivot (2026-07-18).** This document supersedes the earlier
> namelist/JSON-driven `RegcoilProblem` / `RegcoilResult` design. The package is
> now a **scriptable Python object model** with a **stateless** Fortran extension
> holding only the two genuinely expensive kernels. See ADR-019, ADR-020,
> ADR-021, ADR-022. Completed phases 4–5 wrapped the old in-Fortran solve behind
> an opaque handle; that boundary is being re-slimmed in phases 6+ (see
> [PHASES.md](PHASES.md)).

## Package name

Distribution and import name: **`regcoil`** (ADR-015).

## Design constraints (from OVERVIEW + attached design session)

- **No input-file format.** Users drive the calculation with Python scripts and
  set parameters in code. No namelist reader, no JSON schema on the supported
  path (ADR-019). `f90nml` is dropped from the dependency list.
- **Object model, not a monolithic driver.** Three user-facing objects:
  `PlasmaSurface`, `CoilSurface`, `Regcoil`; plus a `Solution` result and a
  `Surface` / `FourierSurface` base (see below). Existing initialization methods
  (VMEC, nescin, FOCUS, circular torus, uniform offset, …) survive as **alternate
  constructors** (classmethods), and the `Surface` ABC keeps this extensible.
- **Stateless Fortran boundary.** `regcoil_variables.f90` is deleted. Fortran is
  2–3 **pure** subroutines with `intent(in)`/`intent(out)` arrays whose extents
  are passed explicitly; no module globals, no derived types crossing the
  boundary, no persistent handle, **no `stop`** (return `info`), and **no
  LAPACK** (ADR-020).
- **All linear algebra in numpy/scipy.** Surface evaluation, basis functions,
  matrix products, the regularized solve, and the λ scan are numpy/scipy. numpy's
  BLAS/LAPACK *is* the same LAPACK, already threaded; keeping it in Python shrinks
  the Fortran layer to something unit-testable (ADR-021).
- **Free λ family.** With `matrix_B`, `matrix_K` symmetric (`matrix_K` positive
  definite), one generalized eigendecomposition diagonalizes the whole λ family;
  scans and target searches become closed-form / cheap. Replaces the Fortran
  auto-regularization + Brent (`regcoil_auto_regularization_solve`, `regcoil_fzero`,
  `regcoil_lambda_scan`).
- **Resolved open decisions (from the design session):**
  1. **Keep** the non-stellarator-symmetry option (boolean `stellarator_symmetric`).
  2. **Remove** Laplace–Beltrami regularization → second derivatives / `nderiv`
     drop out of surface evaluation and matrix assembly (ADR-022).
  3. **Save/load is Python-side NetCDF-4 via `h5netcdf`** (object serialization;
     ADR-028) — see [Saving and loading](#saving-and-loading).
  4. **`Regcoil` accepts a `Surface`** (documented arrays it reads let an outer
     optimization loop pass its own geometry).
- **Binding / build:** `iso_c_binding` + thin C extension `regcoil._core`
  (ADR-017), built with `meson-python` (ADR-002).
- **Deps:** numpy, scipy, matplotlib, plotly; `h5netcdf` for save/load (ADR-028);
  `scipy.io`/optional `netCDF4` for VMEC `wout` reading (ADR-004); pytest for
  tests. Fortran needs a compiler, OpenMP, and BLAS (for `build_g_and_h`'s
  DGEMM) -- no LAPACK, no NetCDF. No SIMSOPT/DESC.

## What a driver script looks like

```python
import numpy as np
from regcoil import PlasmaSurface, CoilSurface, Regcoil

plasma = PlasmaSurface.from_wout("wout_w7x.nc", ntheta=64, nzeta=64, mesh="half")
plasma.set_bnormal_from_bnorm_file("bnorm.w7x")     # or: plasma.Bnormal = <array>

coil = CoilSurface.from_uniform_offset(plasma, separation=0.5,
                                       ntheta=64, nzeta=64, mpol=24, ntor=24)
# alternatives, same downstream:
#   CoilSurface.from_nescin("nescin.w7x", ntheta=64, nzeta=64)
#   CoilSurface.circular_torus(R0=5.5, a=1.5, nfp=5, ntheta=64, nzeta=64)

prob = Regcoil(
    plasma, coil,
    mpol_potential=12, ntor_potential=12,
)                                        # <- all Fortran work happens here

sol = prob.solve(lam=1e-16)
print(sol.f_B, sol.f_K, sol.max_K, sol.rms_K)

scan = prob.scan(np.logspace(-19, -13, 40))   # a SolutionScan (free after the eigendecomposition)
sol  = prob.solve_for_target("max_K", 8.0e6)  # bisection on the same object

K = sol.current_density()          # (3, ntheta_coil, nzeta_coil)
Phi = sol.current_potential()      # (ntheta_coil, nzeta_coil)

# Save the whole run (problem + its plasma & coil + every solution) in one line:
import regcoil
regcoil.save("regcoil_out.w7x.nc", solutions=scan)

# ... later, in a plotting script that never touches the Fortran kernel:
data = regcoil.load("regcoil_out.w7x.nc")
data.plasma.plot()                              # geometry recomputed from Fourier modes
plt.plot(data.solutions.f_K, data.solutions.f_B)  # Pareto front, straight from stored scalars
data.solutions[7].current_potential()           # stored grid, no recompute
```

Three objects, one expensive constructor, cheap solves.

**Immutability.** `Regcoil.__init__` does all the assembly, and the object is
thereafter immutable. Mutating `coil.rmnc` does **not** update a live `Regcoil`;
build a new one. No lazy caches with invalidation — the expensive step is exactly
the step that depends on the geometry, so there is nothing to gain from caching it.

## Python object model

### `Surface` (ABC)

The contract is *evaluation*, not representation:

```python
class Surface(ABC):
    nfp: int
    stellarator_symmetric: bool
    ntheta: int; nzeta: int            # nzetal = nzeta * nfp
    standard_toroidal_angle: bool      # zeta is atan2(y, x)? (ADR-025)

    @abstractmethod
    def _evaluate(self, theta, zetal) -> dict:
        """Return 'r', 'drdtheta', 'drdzeta', each of shape
        (3, ntheta, nzetal), in Cartesian components."""
```

The base class supplies everything derived, identically for all subclasses:
`normal` (cross product), `norm_normal`, `area`, `volume`, `dtheta`, `dzeta`,
the `(theta, zeta, zetal)` grids, and plotting. This is the extensibility hook —
a spline surface, Bezier patch, or coordinate-transformed surface subclasses
`Surface`, implements `_evaluate`, and every downstream consumer works unchanged.

`standard_toroidal_angle` is `True` iff a constant-`zeta` slice of `r` is a
constant *physical* toroidal-angle (`atan2(y, x)`) cross section. Every
constructor below sets/defaults it `True` except
`CoilSurface.from_uniform_offset(..., standard_toroidal_angle=False)` (the
default there — ADR-025): moving plasma points along their normal without
re-solving for the standard toroidal angle leaves `atan2(y, x) != zeta` in
general, so a future cross-section plot (Phase 10) must check this attribute
before assuming constant-`zeta` == constant-physical-angle. Such a surface is
still represented in Fourier space via a **toroidal-angle-shift** mode set
`nu` (`phi = zeta + nu`; ADR-026) — see `FourierSurface`.

> Second derivatives are gone with Laplace–Beltrami regularization (ADR-022), so
> `_evaluate` has no `nderiv` and returns first derivatives only.

### `FourierSurface(Surface)`

The concrete workhorse. Holds `mnmax, xm, xn, rmnc, rmns, zmnc, zmns` and the
optional toroidal-angle-shift modes `numns, numnc` (default zero); its
`_evaluate` is the numpy gemm (`rmnc @ cos(angle)` etc.), placing the point at
`phi = zeta + nu` (`nu` defaults to 0 → `phi = zeta`; ADR-026). The old
`geometry_option_*` integer codes become alternate constructors:

| Legacy option | New constructor |
|---|---|
| plasma 0,1 / coil 0,1 | `circular_torus(R0, a, nfp, ...)` |
| plasma 2,3 | `from_wout(wout, mesh="full"\|"half")` |
| plasma 4 | `from_wout(wout, straight_field_line=True)` |
| plasma 5 | (removed, do not implement in python) |
| plasma 6 | `from_ascii_table(file)` |
| plasma 7 | `from_focus(file)` (also returns Bnormal modes) |
| coil 2, 4 | `from_uniform_offset(plasma, separation, standard_toroidal_angle=False, ...)` |
| coil 3 | `from_nescin(file)` |

`from_uniform_offset`'s `standard_toroidal_angle` argument (default `False`,
ADR-025) selects the construction; it does not distinguish legacy `coil 2` vs
`coil 4` -- the constant-arclength theta reparametrization that separates them
is not implemented (ADR-024 item 4), so both map to the same constructor call
(`standard_toroidal_angle=True` reproduces the `coil 2`/`coil 4` root-solve
behavior for regression purposes). A future `uniform_arclength` argument to
add the `coil 4` reparametrization is a documented idea, not yet implemented.

### `PlasmaSurface(FourierSurface)`

Adds physics attached to the plasma boundary:

- `Bnormal_from_plasma_current` — set from a bnorm file
  (`set_bnormal_from_bnorm_file`), FOCUS modes, or a user array.
- `net_poloidal_current_Amperes`, `curpol`.

### `CoilSurface(FourierSurface)`

Nearly bare — mainly carries `from_uniform_offset` and coil-side Fourier
filtering (`mpol_coil_filter` / `ntor_coil_filter`). `from_uniform_offset`
only calls Fortran when `standard_toroidal_angle=True`; its default
(`False`) is pure Python/numpy (ADR-025).

### `Regcoil`

Constructed from a `plasma`, a `coil`, and solver parameters
(`mpol_potential`, `ntor_potential`, `net_poloidal_current`,
`net_toroidal_current`, `stellarator_symmetric`). Holds the potential Fourier modes
(`xm_potential`, `xn_potential`, `nbf`), `basis_functions`, `g`, `h`,
`matrix_B`, `matrix_K`, `RHS_B`, `RHS_K`, and the cached generalized
eigendecomposition. Methods:

| Method | Role |
|--------|------|
| `solve(lam)` | one regularized solve → `Solution` |
| `scan(lambdas)` | vectorized over a λ array (free after eigendecomposition) |
| `solve_for_target(metric, value)` | bisection/Newton on the closed-form metric |

### `Solution`

A frozen dataclass, one per λ: `lam`, `solution` (mode amplitudes),
`single_valued_current_potential_mn`, `f_B`, `f_K`, `max_K`, `rms_K`,
`max_Bnormal`, `Bnormal_total`, plus **lazy** `current_potential()` /
`current_density()` that expand to grids on demand (keeping grid-sized arrays
lazy matters when scanning many λ). `save(path)` writes the object to
NetCDF-4 (see [Saving and loading](#saving-and-loading)).

A **λ scan** returns a `SolutionScan`, a `Sequence[Solution]` that iterates and
indexes as ordinary `Solution` objects **and** exposes the per-λ scalars as
columnar arrays (`.lam`, `.f_B`, `.f_K`, `.max_K`, `.rms_K`, `.max_Bnormal`)
for direct Pareto/scan plotting. `Regcoil.scan()` and `regcoil.load()` both
return it (ADR-028).

## The solve, in numpy/scipy

`matrix_B` and `matrix_K` are symmetric with `matrix_K` positive definite, so one
generalized eigendecomposition (done once in `Regcoil.__init__`)

```python
w, V = scipy.linalg.eigh(matrix_B, matrix_K)   # ~nbf x nbf (e.g. 600), ~0.1 s
```

diagonalizes the whole family:
`x(λ) = V @ ((Vᵀ·b_B + λ·Vᵀ·b_K) / (w + λ))`. Every λ afterwards is O(nbf²), and
`f_B(λ)`, `f_K(λ)` are smooth closed-form scalar functions to bisect or
Newton on directly — replacing `regcoil_auto_regularization_solve.f90`,
`regcoil_fzero.f`, and `regcoil_lambda_scan.f90` with ~15 lines of Python.

## Saving and loading

All I/O is Python (ADR-028). Files are **NetCDF-4 written via `h5netcdf`** — HDF5
underneath (no NetCDF C library in the wheels/CI), but with named, shared
dimensions so the file is self-documenting and directly `xarray.open_dataset`-able.
There is no attempt to preserve the legacy Fortran `regcoil_out.*.nc` layout; this
is a fresh Python-defined schema.

### API

One polymorphic module-level pair, plus thin instance methods:

```python
import regcoil

# One line saves the whole run. `solutions` transitively carries its
# problem -> coil, plasma, so all four object kinds land in one file:
regcoil.save("regcoil_out.w7x.nc", solutions=scan)
regcoil.save("regcoil_out.w7x.nc", problem=prob, solutions=scan)  # explicit, same file

# Any subset is fine; only the groups you pass are written:
plasma.save("plasma.nc")          # == regcoil.save("plasma.nc", plasma=plasma)
prob.save("prob.nc"); sol.save("sol.nc")

# load() returns a container; any attribute may be None if that group is absent:
data = regcoil.load("regcoil_out.w7x.nc")
data.plasma, data.coil, data.problem, data.solutions   # solutions is a SolutionScan

# Per-class load returns the single object (errors if that group is missing):
PlasmaSurface.load("plasma.nc");  Regcoil.load("prob.nc")
```

`regcoil.save(path, *, plasma=None, coil=None, problem=None, solutions=None)`
writes only the groups it is given and **dedupes the shared problem/surfaces** —
an N-λ scan stores its `Regcoil` (and plasma, coil) once, not N times.

### File layout

Groups by object, each tagged with a `_class` attribute; the root carries a
`format_version` for dispatch and future migration:

```
/                format_version, regcoil __version__
/plasma          _class=PlasmaSurface : xm,xn,rmnc,rmns,zmnc,zmns,numns,numnc,
                   nfp,ntheta,nzeta,stellarator_symmetric,standard_toroidal_angle,
                   net_poloidal_current, curpol, Bnormal_from_plasma_current
/coil            _class=CoilSurface   : the FourierSurface fields (incl. nu modes)
/problem         mpol_potential, ntor_potential, net_poloidal_current,
                   net_toroidal_current, stellarator_symmetric,
                   Bnormal_from_net_coil_currents(ntheta_plasma,nzeta_plasma)
                   # ^ stored (ADR-029): λ-independent, its recompute is the
                   #   Fortran h term, so the "net coil currents" Bnormal panel
                   #   plots from a saved file with no kernel call
/solutions       dim lambda(nlambda):
                   lam, solution(nlambda,nbf),
                   f_B, f_K, max_K, rms_K, max_Bnormal,
                   Bnormal_total(nlambda,ntheta_plasma,nzeta_plasma),
                   current_potential(nlambda,ntheta_coil,nzeta_coil),
                   current_density(nlambda,3,ntheta_coil,nzeta_coil)
```

A λ scan is **flattened along a `lambda` dimension** (a scan is homogeneous by
construction — one shared problem, one `nbf`, one grid — and analysis reads want
columnar arrays); a single `solve()` is the `nlambda=1` case.

### What is stored vs recomputed

The rule (ADR-028): **store the source of truth, plus exactly those derived
quantities whose recomputation would need the Fortran kernel or the big operators
(`g`, `h`, `matrix_B/K`, `V`); recompute everything else from the surfaces and the
mode amplitudes.** The big operators are **never** saved.

- **Surfaces** store only their Fourier source of truth (plus the plasma's
  `Bnormal_from_plasma_current`, `net_poloidal_current`, `curpol`, which are not
  recoverable from the coefficients). `r`, `normal`, `area`, `volume` are one gemm
  away, rebuilt lazily on load.
- **`Regcoil`** stores its scalar parameters and references to plasma/coil.
  `basis_functions`, `f_all`, `d_xyz` are cheap numpy, rebuilt on load; `g`, `h`,
  `matrix_B/K`, `RHS_*`, `w`, `V` are not stored. **One exception (ADR-029):** the
  λ-independent plasma-grid array `Bnormal_from_net_coil_currents` *is* stored (it
  derives from the Fortran `h` term), so the Bnormal "from net coil currents" plot
  panel needs no kernel on a loaded run.
- **`Solution`** stores `lam`, `solution` (amplitudes), the scalar diagnostics,
  **and** the result-sized grids `Bnormal_total`, `current_potential`,
  `current_density`. `Bnormal_total` *must* be stored (its recompute is
  `g @ solution`); the other two are stored because recomputing them would
  materialize the ~0.5 GB `f_all` intermediate at 128×128, whereas the result
  grids are ~0.4 MB/λ. Storing them makes every plot target available with **no
  Fortran call, no expensive BLAS, and no object reconstruction** — the goal a
  saved file exists to serve.

### Loaded objects

`Regcoil.__init__` is split into a **cheap** assembly (surfaces, potential modes,
`basis_functions`, `f_all`, `d_xyz`) and an **expensive** one (`build_g_and_h`,
`matrix_B/K`, `eigh`). `load()` runs only the cheap part, so a loaded problem's
Solutions plot instantly; the operators are rebuilt lazily (with a log message)
only if the user asks for a *new* λ via `solve()`/`scan()` on the loaded object.

## Plotting and visualization

Phase 10; architecture in ADR-029. The design goals: make **one plot type at a
time** as easily as the full dashboard; plot a **live** run or a **saved** run
identically and with **no Fortran kernel / no expensive BLAS**; and let plots be
composed (subplots, overlays) freely.

### Principles

1. **Object model only.** Plot functions take in-memory objects
   (`PlasmaSurface`, `CoilSurface`, `Regcoil`, `Solution`, `SolutionScan`) or the
   `regcoil.load()` container — never file paths or `xarray`. `load()` rebuilds
   real objects (ADR-028 cheap-assembly split), so a live and a loaded object are
   the same interface. File → object is the CLI's job.
2. **Hybrid placement.** Canonical free functions live in `regcoil.plot`; thin
   `obj.plot_*()` methods delegate to them. matplotlib and plotly are imported
   **lazily inside** the functions, so `import regcoil` pulls in no plotting
   library and the core solve keeps no plotting dependency.
3. **Atomic functions return their axes/figure; dashboards compose them.** Every
   atomic function accepts `ax=` (matplotlib) or `fig=` (plotly) and returns it.
   `plot.all()` and the λ-scan grids only build the layout and call the atomic
   functions — no plotting logic of their own.
4. **matplotlib for 2D, Plotly for 3D interactive.**
   `regcoil.plot.DEFAULT_FIGSIZE = (14.5, 8.1)` is used whenever a matplotlib
   function must create its own figure (i.e. no `ax=` was passed).

### Functions

```python
import numpy as np
import regcoil
from regcoil import plot

data = regcoil.load("regcoil_out.w7x.nc")   # or use live objects from a script

# --- 2D, matplotlib ---

# Surface cross sections. Default phi = np.array([0, 0.5, 1, 1.5]) * pi / nfp
# (half a field period); pass any phi array. Uses Surface.cross_section(phi),
# correct for both standard_toroidal_angle values (ADR-025).
plot.cross_section(data.plasma, data.coil)
plot.cross_section(data.plasma, data.coil, phi=np.linspace(0, np.pi/data.plasma.nfp, 6))

# Pareto front — one or several runs overlaid; choose the axes.
plot.pareto(data.solutions, x="f_K", y="f_B")
plot.pareto([scanA, scanB], x="max_K", y="max_Bnormal", labels=["A", "B"])

# lambda traces (f_B, f_K vs lambda), the rest of the legacy scan grid.
plot.lambda_scan(data.solutions)

# (theta, zeta) field maps — one lambda (Solution) at a time…
plot.current_potential(data.solutions[7], kind="single_valued")   # or "total"
plot.current_density(data.solutions[7])
plot.bnormal(data.solutions[7], component="total")   # "plasma_current"|"net_coil"|"total"
# …or the multi-lambda panel grid over a scan (composition of the atomic fn):
plot.current_potential_scan(data.solutions, kind="total", nmax=16)
plot.bnormal_scan(data.solutions, nmax=16)

# --- 3D interactive, Plotly ---

# One scene, ANY subset of {plasma, winding surface, cut coils}; omit what you
# don't want. The winding surface renders as "wireframe", "translucent", or "solid".
plot.plot_3d(plasma=data.plasma, winding_surface=data.coil, winding_surface_style="wireframe")
plot.plot_3d(winding_surface=data.coil, coils=cut)        # winding surface + cut coils
plot.plot_3d(plasma=data.plasma, winding_surface=data.coil, coils=cut,
             winding_surface_style="translucent")

# --- everything (pure composition of the above) ---
plot.all(data)                       # ≈ the legacy regcoilPlot dashboard
```

Every function also has an object-method alias, e.g. `data.plasma.plot_cross_section(...)`,
`data.solutions.plot_pareto(...)`, `data.solutions[7].plot_current_potential(...)`.

### Coil cutting (compute, separate from plotting)

Cutting is a computation that produces coil geometry; plotting it is downstream
(ADR-029). The irreversible MAKEGRID file write stays out of the plot path.

```python
cut = regcoil.cut(data.solutions[7], coils_per_half_period=3, thickness=0.05)
cut.save_makegrid("coils.w7x")       # replaces cutCoilsFromRegcoil's file output
plot.plot_3d(winding_surface=data.coil, coils=cut)     # or plot.coil_3d(cut)
```

`cut()` contours the total current potential, maps contours to 3D curves via the
coil Fourier modes, and (given `thickness`) builds the finite-thickness ribbons —
the Plotly port of `m20160811_01_plotCoilsFromRegcoil.m`.

## Command-line interface

One console entry point, `regcoil`, with subcommands (ADR-029, ADR-015). Each
subcommand is a thin wrapper: `regcoil.load(file)` → plot/cut → `show`/save. It is
the only layer that touches file paths; the library functions stay object-only.

```text
regcoil plot    FILE.nc [--only NAME]... [--phi ...] [--lambda-index N]
                        [--save OUT.(png|html)] [--no-show]
        # No --only: the full dashboard (plot.all). --only selects atomic plots,
        # repeatable: cross_section | pareto | lambda_scan | current_potential |
        # current_density | bnormal | plot_3d.  --phi overrides cross-section
        # angles; --lambda-index picks which Solution for the (theta,zeta) maps.

regcoil compare FILE.nc [FILE.nc ...] [--x f_K|max_K] [--y f_B|max_Bnormal]
                        [--save OUT.png] [--no-show]
        # Overlays plot.pareto across several saved runs (replaces compareRegcoil).

regcoil cut     FILE.nc --coils-per-half-period N [--lambda-index N]
                        [--thickness T] [--out coils.NAME] [--plot] [--no-show]
        # Runs regcoil.cut(...); writes a MAKEGRID coils.* file; --plot shows the
        # 3D finite-thickness figure (replaces cutCoilsFromRegcoil / cut_saddle_coil).
```

Notes:
- 3D subcommands (`plot_3d`, `cut --plot`) open an interactive Plotly view, or
  write a standalone `.html` with `--save`. 2D `--save` writes an image.
- `--no-show` (for headless/CI/scripted image dumps) skips the interactive window.
- No `regcoilPlot` / `compareRegcoil` console names — one `regcoil` CLI only.

## Fortran extension boundary

The whole boundary is 2 mandatory functions plus 1 optional — no shared state,
no derived types crossing it.

### Stays in Fortran

- **(a) Pairwise kernel** — the only genuinely O(N²) hot path.
- **(b) Uniform-offset surface** — serial root solves + splines + a DFT.
- **(c, optional) Surface evaluation** — a gemm; start in numpy, add here only if
  profiling demands it.

### Moves to Python (numpy/scipy)

- Surface evaluation from Fourier coefficients, `normal`, `area`, `volume`, grids.
- `basis_functions`, `f_x/f_y/f_z`, `d_x/d_y/d_z`.
- `matrix_B = gᵀ(g/N)`, `matrix_K = Σ fᵢᵀ(fᵢ/N)`, the RHS vectors.
- The regularized solve, the generalized eigendecomposition, the λ scan, and any
  target/auto-λ search.
- All I/O: VMEC `wout`, nescin, bnorm, FOCUS, and NetCDF output.
- Plot / compare / cut / Plotly tools; string-option normalization.

### Deleted (not ported)

- `regcoil_variables.f90` (module globals) and the opaque-handle C API built in
  Phase 5.
- The in-Fortran solve chain: `regcoil_prepare_solve.f90`, `regcoil_solve.f90`,
  `regcoil_diagnostics.f90`, and the globals-coupled matrix build.
- `regcoil_auto_regularization_solve.f90`, `regcoil_lambda_scan.f90`,
  `regcoil_compute_lambda.f90`, `regcoil_fzero.f` (λ search moves to scipy;
  audit whether any geometry root solve is still needed by the offset routine).
- Adjoint / sensitivity / WSO; `regcoil_svd_scan`; Laplace–Beltrami metric block.
- Fortran namelist read/write and NetCDF (`ezcdf`, `read_wout` NetCDF path).
- EFIT input format.

## The interface signatures

Presented as pure Fortran subroutines (the mathematical contract). Under
`iso_c_binding` (ADR-017), a thin `bind(C)` wrapper in `fortran/regcoil_c_api.f90`
passes the extents explicitly and marshals the NumPy buffers; `src/regcoil/_core.c`
exposes them to Python and turns a nonzero `info` into an exception.

### (a) Pairwise kernel — `regcoil_build_g_and_h`

Fuse `inductance @ basis_functions` inside Fortran so the 134 MB inductance matrix
is never materialized:

```fortran
subroutine regcoil_build_g_and_h( &
     ntheta_plasma, nzeta_plasma, ntheta_coil, nzeta_coil, nfp, nbf, &
     r_plasma, normal_plasma, r_coil, normal_coil, drdtheta_coil, drdzeta_coil, &
     basis_functions, net_poloidal_current, net_toroidal_current, &
     dtheta_coil, dzeta_coil, g, h, info)

  integer,  intent(in)  :: ntheta_plasma, nzeta_plasma, ntheta_coil, nzeta_coil, nfp, nbf
  real(dp), intent(in)  :: r_plasma(3, ntheta_plasma, nzeta_plasma)
  real(dp), intent(in)  :: normal_plasma(3, ntheta_plasma, nzeta_plasma)
  real(dp), intent(in)  :: r_coil(3, ntheta_coil, nzeta_coil*nfp)
  real(dp), intent(in)  :: normal_coil(3, ntheta_coil, nzeta_coil*nfp)
  real(dp), intent(in)  :: drdtheta_coil(3, ntheta_coil, nzeta_coil*nfp)
  real(dp), intent(in)  :: drdzeta_coil(3, ntheta_coil, nzeta_coil*nfp)
  real(dp), intent(in)  :: basis_functions(ntheta_coil*nzeta_coil, nbf)
  real(dp), intent(in)  :: net_poloidal_current, net_toroidal_current
  real(dp), intent(in)  :: dtheta_coil, dzeta_coil
  real(dp), intent(out) :: g(ntheta_plasma*nzeta_plasma, nbf)
  real(dp), intent(out) :: h(ntheta_plasma*nzeta_plasma)
  integer,  intent(out) :: info
end subroutine
```

Blocked internally: accumulate `inductance` for a chunk of plasma rows into a
small buffer, DGEMM that chunk against `basis_functions` (faster than the
`matmul` intrinsic), write into `g`. Keeps the working set in cache and drops
peak memory ~100×. Add OpenMP over plasma rows and mark the wrapper
threadsafe so the GIL is released during the kernel.

### (a′) Debug variant — `regcoil_build_inductance`

Same argument list minus `basis_functions`, returning the full `inductance`
matrix, as a **separate** entry point for debugging and regression tests against
the legacy code. Two subroutines is cleaner than one with `optional` outputs.

### (b) Offset surface — `regcoil_uniform_offset_surface`

Crucially returns **Fourier coefficients**, so a uniform-offset coil is just a
`FourierSurface` like any other and nothing downstream special-cases it:

```fortran
subroutine regcoil_uniform_offset_surface( &
     mnmax_in, xm_in, xn_in, rmnc_in, rmns_in, zmnc_in, zmns_in, lasym, nfp, &
     separation, mpol_out, ntor_out, ntheta_transform, nzeta_transform, tol, &
     mnmax_out, xm_out, xn_out, rmnc_out, rmns_out, zmnc_out, zmns_out, info)
```

`mnmax_out = mpol_out*(2*ntor_out+1) + ntor_out + 1` is known from the inputs, so
the output arrays are sized deterministically — no allocatable-return gymnastics.
Internally: sample a fine grid on the offset surface, root-solve back to the
plasma surface per point (splines for interpolation), then DFT to coefficients.

### (c) Optional — `regcoil_evaluate_surface`

`regcoil_evaluate_surface(mnmax, xm, xn, rmnc, ..., theta, zetal, r, drdtheta,
drdzeta, ..., info)`. Start with numpy; add this only if profiling says so.

## Conventions to pin down (assert in the Python wrapper)

These are where a rewrite bleeds time — document them at the top of the Fortran
module and assert them at the boundary:

- **Array order.** Fortran is column-major; numpy defaults to row-major. Allocate
  boundary arrays Python-side with `order="F"` (or `np.asfortranarray`) and assert
  `arr.flags.f_contiguous` so a copy is a loud choice, not a hidden multi-MB cost.
- **Index order** is `(3, ntheta, nzetal)` with the Cartesian component first
  (matches the legacy code and keeps the innermost loop contiguous).
- **`zeta` vs `zetal`.** Plasma quantities live on one field period (`nzeta`); coil
  quantities span the full torus (`nzetal = nzeta*nfp`). Encode it in argument
  names and check extents.
- **`xn` includes the factor of `nfp`** (VMEC convention) and the angle is
  `m·θ − n·ζ`. FOCUS files omit it, so apply `xn = xn * nfp` **inside** the
  `from_focus` constructor, keeping the invariant everywhere else.
- **Units:** SI, with `mu0/(4π)` folded into the inductance and
  `dtheta_coil*dzeta_coil*mu0/(8π²)` into `h`, as in the legacy code.

## Build system

`meson-python` backend (ADR-002); `iso_c_binding` + `regcoil._core` C extension
(ADR-017). The extension links the Fortran runtime, OpenMP, and BLAS (DGEMM in
`build_g_and_h`) — **no LAPACK, no NetCDF** — after the boundary is
re-slimmed. `numpy.distutils` is gone as of Python 3.12, so no f2py path.

## Testing

- **pytest** for all Python tests and for driving the Fortran kernels via
  `regcoil._core` with small manufactured inputs (ADR-008).
- Layout: `tests/unit/`, `tests/integration/`; examples wired as regression tests
  against the physics tolerances (the contract, not binary-identical NetCDF).
- Key numeric checks: `build_g_and_h` vs `build_inductance @ basis_functions`;
  the offset-surface coefficients vs the legacy routine; the closed-form λ family
  vs direct per-λ solves; save/load round-trip against the ADR-028 schema.
