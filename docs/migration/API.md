# Target API sketch

Design target for the migration—not a frozen public API. Refine via [DECISIONS.md](DECISIONS.md).

## Package name

Working name: **`regcoil`** (see ADR-015).

## Design constraints (from OVERVIEW)

- Multiple in-memory instances; **no Fortran globals** for problem state.
- Qualitative options are **strings**.
- Input: **namelist via `f90nml`** and **JSON**.
- No SVD scan; no adjoint/WSO.
- Tools (plot, compare, cut coils) live in the package.
- Deps: numpy, scipy, matplotlib, f90nml, plotly, optional netCDF4; pytest for tests; BLAS/LAPACK for Fortran. No SIMSOPT/DESC.

## Python surface

### High-level

```python
from regcoil import run, RegcoilProblem, RegcoilResult

result = run("regcoil_in.my_case")          # namelist (f90nml)
result = run("my_case.json")                # JSON
problem = RegcoilProblem.from_file("my_case.json")
# or RegcoilProblem.from_namelist(...) / from_dict(...)
result = problem.solve()
result.write_netcdf("regcoil_out.my_case.nc")
```

Two problems must not share mutable Fortran state:

```python
a = RegcoilProblem.from_file("a.json")
b = RegcoilProblem.from_file("b.json")
ra, rb = a.solve(), b.solve()  # independent
```

### String options (qualitative)

Integers that only select among qualitatively different modes become strings. Counts, resolutions, and indices stay numeric.

Illustrative mapping (exact spellings TBD in docs / ADR-009; align namelist, JSON, and Python):

| Legacy field | Legacy values | Proposed string values (examples) |
|--------------|---------------|-------------------------------------|
| `general_option` | 1, 2, 4, 5 (3 **removed**) | `"lambda_scan"`, `"nescout_potential"`, `"auto_regularization"`, `"auto_regularization_bounded"` (names TBD) |
| `geometry_option_plasma` | 0, 1, 2, … | `"circular_torus"`, `"vmec"`, … |
| `geometry_option_coil` | 0, 1, … | `"offset"`, `"nescin"`, … |
| `symmetry_option` | 1, 2, … | `"stellarator"`, … |
| `target_option` | already partly strings | keep / normalize string enums |
| `regularization_term_option` | already strings | keep |

Validation should reject unknown strings with a clear error. Policy for old integer codes in namelists: ADR-009.

### `RegcoilProblem`

Owns configuration and a handle to instance Fortran state:

- Resolution, geometry, currents, regularization, verbosity
- Methods:

| Method | Role |
|--------|------|
| `from_file` / `from_namelist` / `from_dict` | Load namelist (`f90nml`) or JSON |
| `to_json` / `to_namelist` | Round-trip helpers for tests/docs |
| `load_geometry()` | Python I/O → arrays |
| `build()` | Fortran matrix assembly + prepare (instance) |
| `solve_lambda(lambda_)` | One regularized solve + diagnostics |
| `scan_lambda(lambdas)` | Former option 1 |
| `find_lambda(...)` | Former options 4/5 via SciPy |
| `solve()` | Dispatch on string run mode |

**Not present:** `svd_scan`.

### `RegcoilResult`

Metrics and fields needed by tests and tools (`lambda`, `chi2_B`, `chi2_K`, `max_K`, `max_Bnormal`, `exit_code`, potentials, grids, …). NetCDF write/read via the library chosen in ADR-004.

### Residual for auto-regularization

```python
def residual(log_lambda: float, problem: RegcoilProblem) -> float:
    """Physics residual only: solve at λ and return (metric - target)."""
    ...
```

SciPy root finder stays outside this function.

### CLI

```bash
regcoil run regcoil_in.my_case
regcoil run my_case.json
regcoil plot regcoil_out.my_case.nc
regcoil compare regcoil_out.a.nc regcoil_out.b.nc
regcoil cut-coils ...
regcoil plot-coils ...   # Plotly port of m20160811_01_*.m
```

Exact subcommand names TBD; console scripts listed in `pyproject.toml`.

### Plotting / coil tools (in-package)

| Legacy | Package role |
|--------|----------------|
| `regcoilPlot` | matplotlib diagnostics plots |
| `compareRegcoil` | multi-run comparison plots |
| `cutCoilsFromRegcoil`, `cut_saddle_coil` | contour → filament coils |
| `m20160811_01_plotCoilsFromRegcoil.m` | Plotly 3D coil visualization |

## Input schema

- **Namelist:** same logical `&regcoil_nml` (or successor) fields; qualitative fields are character strings; parsed with `f90nml` in Python (Fortran namelist reader retired for the supported path).
- **JSON:** object mirroring namelist keys (plus optional nested groups later). Example:

```json
{
  "general_option": "lambda_scan",
  "geometry_option_plasma": "vmec",
  "nlambda": 20,
  "mpol_potential": 8,
  "ntor_potential": 8
}
```

## Fortran extension boundary

### Stays in Fortran

- Matrix build, prepare, solve, per-λ diagnostics
- Instance derived type (or explicit state args)—**no problem globals**
- BLAS/LAPACK

### Moves to Python

- Driver, namelist/JSON I/O, validation orchestration
- λ grid, scans, auto-λ root search
- NetCDF read/write; VMEC NetCDF wout ingest
- Plot / compare / cut / Plotly tools
- String option normalization

### Deleted (not ported)

- Adjoint / sensitivity / WSO
- `regcoil_svd_scan` / option 3
- All MATLAB except the temporary Plotly source

### Transitional gray area

Globals and Fortran namelist read may remain until Phases 5–6. Each must be removed before Phase 12.

## Binding style

Build backend is **meson-python** (ADR-002). Binding style is **`iso_c_binding`** + thin C extension `regcoil._core` (ADR-017). Requirement: NumPy in/out as needed; instance handle or state object in Phase 5; no NetCDF in the solve path after Phase 8.

## Testing

- **pytest** for all Python tests and for driving Fortran-backed unit tests unless ADR-008 chooses a Fortran framework.
- Layout: `tests/unit/`, `tests/integration/`, examples wired as regression tests.
- Concurrent-instance test is mandatory after Phase 5.
