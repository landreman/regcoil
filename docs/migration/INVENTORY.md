# Codebase inventory (pre-migration)

Snapshot to support kill lists and scoping. Prefer updating PHASES checklists over continually editing this file.

## Top-level layout (current)

| Path | Role | Migration fate |
|------|------|----------------|
| `fortran/regcoil.f90` | Fortran `PROGRAM` driver | Retire (Phase 12); logic → Python |
| `fortran/regcoil_*.f90`, `fortran/regcoil_fzero.f` | Core sources | Library under `fortran/`; deglobalize |
| `fortran/regcoil_variables.f90` | **Module globals** | Replace with derived type / args (Phase 5) |
| `fortran/mini_libstell/` | kinds, wout, ezcdf, … | Trim; drop NetCDF/ezcdf |
| `makefile`, `makefile.depend` | Host-specific build | Packaging replaces for users |
| `examples/` | Regressions | Keep; drive with pytest |
| `equilibria/` | Sample inputs | Keep |
| `manual/` | LaTeX manual | Replace with RTD (Phase 11) |
| `windingSurfaceOptimization/` | Surface optimization | **Delete** (Phase 1) |
| `regcoilPlot`, `compareRegcoil`, `cutCoilsFromRegcoil`, `cut_saddle_coil` | Python CLIs | **Fold into package** (Phase 9) |
| `*.m`, `coilMetricScripts/` | MATLAB | **Delete**; port only `m20160811_01_plotCoilsFromRegcoil.m` → Plotly |
| `.github/workflows/publish_manual.yml` | LaTeX → gh-pages | Replace with RTD / pytest CI |

## Fortran module roles

### Driver / orchestration

- `regcoil.f90` — main program
- `regcoil_read_input.f90`, `regcoil_validate_input.f90`, `regcoil_write_input.f90` → Python (`f90nml` + JSON)
- `regcoil_compute_lambda.f90`, `regcoil_lambda_scan.f90` → Python
- `regcoil_svd_scan.f90` — **delete** (ADR-010)
- `regcoil_auto_regularization_solve.f90` — algorithm → Python/SciPy; solve/metrics stay callable
- `regcoil_fzero.f` — audit (auto-λ vs geometry roots)
- `regcoil_compute_diagnostics_for_nescout_potential.f90` — keep behavior via Python driver

### Geometry / setup

- `regcoil_init_plasma_mod.f90`, coil surface / offset / expand, nescin, bnorm, EFIT, Fourier modes, splines
- VMEC wout NetCDF path → Python readers (Phase 8)

### Matrices / linear algebra (keep in Fortran; instance state)

- `regcoil_build_matrices.f90`, `regcoil_prepare_solve.f90`, `regcoil_solve.f90`, `regcoil_diagnostics.f90`
- Today coupled to `regcoil_variables` globals — **must become instance-based**

### I/O

- `regcoil_write_output.f90` + `ezcdf*` → Python NetCDF (ADR-004)
- `mini_libstell/read_wout_mod.F` — NetCDF portions → Python

### Adjoint / sensitivity (**delete**)

- `regcoil_init_sensitivity.f90`, `regcoil_adjoint_solve.f90`, `regcoil_fixed_norm_sensitivity.f90`
- Related variables, validate branches, write_output fields, Fourier sensitivity init
- `windingSurfaceOptimization/`, `manual/adjoint.tex`

### SVD (**delete**)

- `regcoil_svd_scan.f90`
- `general_option == 3` in driver / validate / output / docs

## Qualitative integer options (→ strings)

From `regcoil_variables.f90` (non-exhaustive):

| Variable | Type today | Notes |
|----------|------------|-------|
| `general_option` | integer | 3 = SVD → remove; others → strings |
| `geometry_option_plasma` | integer | → strings |
| `geometry_option_coil` | integer | → strings |
| `symmetry_option` | integer | → strings |
| `sensitivity_option` / `sensitivity_symmetry_option` | integer | delete with adjoint |
| `regularization_term_option` | character | already string-like — keep |
| `target_option` | character | already string-like — keep |

## Globals problem

Almost all subroutines `use regcoil_variables`. Phase 5 replaces this with a derived type (or explicit state) passed through the API so two Python `RegcoilProblem` objects do not clobber each other.

## NetCDF touch points

| Location | Use |
|----------|-----|
| `makefile` host stanzas | `-lnetcdf -lnetcdff` / `nc-config` |
| `regcoil_write_output.f90` | Writes `regcoil_out*.nc` |
| `mini_libstell/ezcdf*`, `read_wout_mod.F` | Fortran NetCDF |
| `examples/testsCommon.py`, `regcoilPlot`, `compareRegcoil`, cut scripts | `scipy.io.netcdf_file` today |

## Lambda search

| Piece | Fate |
|-------|------|
| Staging / Brent in `regcoil_auto_regularization_solve.f90` | Python + SciPy |
| Per-λ solve + metrics | Fortran, instance methods |
| Fixed λ list / scan | Python |

## MATLAB inventory (**delete** except temporary Plotly source)

- `m20160811_01_plotCoilsFromRegcoil.m` — **port to Plotly**, then delete
- `m20160811_02_compare2CoilsetsFromRegcoil.m` — delete
- `regcoil.m` — delete
- `coilMetricScripts/*.m` — delete (no Python port planned)

## Existing tests

- Example discovery via `tests.py` + `regcoil_in.<dirname>`
- Runner: `examples/runExamples.py` under `make test`
- Assertions via `testsCommon` + NetCDF
- Migration: pytest; add **unit** tests (Python + Fortran cores per ADR-008); keep example tolerances as physics contract
- Concurrent-instance test required after deglobalization

## Kill / keep / move summary

| Keep (Fortran core) | Delete | Move to Python package |
|---------------------|--------|------------------------|
| build_matrices, prepare_solve, solve, diagnostics (instanced) | adjoint, sensitivity, WSO | driver, λ search, NetCDF I/O |
| geometry numerics as needed | SVD scan | namelist (`f90nml`) + JSON |
| trimmed mini_libstell (no ezcdf) | all MATLAB after Plotly port | plot, compare, cut-coils, Plotly coils |
| | LaTeX-as-canonical manual | RTD docs |
