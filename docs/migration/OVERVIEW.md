# pyREGCOIL migration overview

Convert REGCOIL from a standalone Fortran executable into a **pip-installable Python package** built around a **scriptable object model** (`PlasmaSurface`, `CoilSurface`, `Regcoil`, `Solution`). Python/numpy own I/O, the surface geometry, **all** matrix products, the regularized solve, the λ scan, plotting/coil tools, and documentation. Fortran shrinks to a **stateless** extension holding only the two genuinely expensive kernels — the O(N²) pairwise inductance/field sum and the uniform-offset surface — as pure `intent(in)`/`intent(out)` functions with **no module globals, no derived types crossing the boundary, no persistent handle, and no LAPACK**.

> **Architecture pivot (2026-07-18).** This overview was revised to match the
> object-model / stateless-kernel design (ADR-019 – ADR-022). It supersedes the
> earlier "namelist/JSON driver over an instance-based Fortran solve" framing.
> See [API.md](API.md) and [PHASES.md](PHASES.md).

Companion docs:

- [PHASES.md](PHASES.md) — ordered work packages and exit criteria
- [API.md](API.md) — intended Python / Fortran boundary
- [DECISIONS.md](DECISIONS.md) — append-only design log
- [INVENTORY.md](INVENTORY.md) — current-code map (roles, kill lists, touch points)
- [LOCAL.md](LOCAL.md) — how to build/run `make test` locally (esp. macOS Homebrew)

## Goals

1. **Python packaging** — `pip install` builds the Fortran extension and exposes the Python API.
2. **Python-driven execution** — remove the Fortran `PROGRAM regcoil` as the user-facing entry point.
3. **Fortran only for the two hot paths** — the O(N²) pairwise inductance/field kernel and the uniform-offset surface stay compiled; everything else is numpy/scipy. (ADR-020)
4. **Stateless Fortran** — delete `regcoil_variables`; entry points are pure `intent(in)`/`intent(out)` functions returning `info` (no `stop`), with no globals, no derived types crossing the boundary, and no persistent handle. (ADR-020)
5. **Linear algebra in numpy/scipy** — matrix products and the regularized solve are numpy/scipy; one generalized eigendecomposition makes the whole λ family closed-form (no Fortran BLAS/LAPACK *for the solve*, no Brent; the unrelated DGEMM call inside `regcoil_build_g_and_h`, ADR-020, is a performance choice for that one kernel, not a solve). (ADR-021)
6. **Scriptable object model** — `PlasmaSurface`, `CoilSurface`, `Regcoil`, `Solution`; existing initialization methods survive as alternate constructors; no input-file format. (ADR-019)
7. **Boolean symmetry switch for `Regcoil`** — `stellarator_symmetric=True` uses the stellarator-symmetric basis; `False` includes both sine and cosine basis families (ADR-019).
8. **NetCDF only from Python** — no NetCDF compile/link dependency in Fortran.
9. **Fold tooling into the package** — `regcoilPlot`, `compareRegcoil`, and coil-cutting scripts become `regcoil` modules / CLI entry points.
10. **Plotly coil plot** — port `m20160811_01_plotCoilsFromRegcoil.m` to Python with Plotly; **remove all other MATLAB** sources.
11. **Remove SVD scan** — delete `regcoil_svd_scan` and `general_option == 3` pathways.
12. **Remove adjoint / winding-surface optimization** — delete sensitivity / adjoint Fortran and `windingSurfaceOptimization/`.
13. **Unit + regression tests** — pytest for Python and integration; unit tests for Fortran cores (framework TBD; see ADR-008). CI on GitHub Actions.
14. **Modern docs** — replace the LaTeX manual with a Read the Docs site (Sphinx or equivalent); **remove** `.github/workflows/publish_manual.yml`.
15. **Layout** — `src/`, `tests/`, `docs/`, `fortran/` (or similar) + `pyproject.toml`.

## Design principles

- **Minimal dependencies.**
  - **Python (allowed):** `numpy`, `scipy`, `matplotlib`, `plotly`, and a NetCDF library (see ADR-004). **Testing:** `pytest`. (`f90nml` dropped with the input-file format — ADR-019.)
  - **Fortran:** a compiler, OpenMP, and BLAS (DGEMM in `regcoil_build_g_and_h`, chosen over the `matmul` intrinsic for performance — ADR-020). **No LAPACK** (numpy/scipy provide it for the solve — ADR-021) and no NetCDF Fortran libs after Phase 9.
  - **Do not** depend on **SIMSOPT**, **DESC**, or other large stellarator stacks.
- **Pytest** is the Python test runner (unit, integration, examples).
- Prefer **instance-friendly** APIs (classes / derived types) over process-global state.
- Prefer explicit typed options for qualitative choices; keep integers for counts, indices, and resolutions.
- One phase per PR when practical; do not mix directory moves with physics refactors.
- Record non-obvious choices in [DECISIONS.md](DECISIONS.md).

## Non-goals (this migration)

- Rewriting the O(N²) pairwise inductance/field kernel or the offset-surface root-solve in Python (they stay in Fortran).
- Preserving the Fortran namelist input file format — parameters are set in Python code (ADR-019).
- Preserving binary-identical NetCDF files if metadata/order changes slightly; **example regression tolerances** are the contract.
- Supporting adjoint-based winding-surface design, SVD scan, or Laplace–Beltrami regularization in the new package (ADR-007, ADR-010, ADR-022).
- Porting `coilMetricScripts/*.m`, `regcoil.m`, or `m20160811_02_compare2CoilsetsFromRegcoil.m` (delete with other MATLAB).
- Introducing SIMSOPT/DESC (or similar) as dependencies.

## Target architecture

```text
pip package: regcoil
│
├── Python (numpy / scipy)
│   ├── object model: Surface → FourierSurface → PlasmaSurface / CoilSurface
│   ├── Regcoil (assembles matrices, one generalized eigendecomposition), Solution
│   ├── geometry / VMEC / nescin / bnorm / FOCUS loaders (string options)
│   ├── basis functions, matrix_B / matrix_K, regularized solve
│   ├── λ family + target search (closed-form after the eigendecomposition)
│   ├── NetCDF I/O (library: ADR-004)
│   ├── plot / compare / cut-coils tools (matplotlib + plotly)
│   └── docs → Read the Docs
│
└── Fortran extension (compiled; runtime + OpenMP + BLAS; no LAPACK, no NetCDF)
      ├── regcoil_build_g_and_h   — fused inductance @ basis → g, h  (pure, info)
      ├── regcoil_build_inductance — full matrix, debug/regression   (pure, info)
      └── regcoil_uniform_offset_surface — returns Fourier coeffs    (pure, info)
```

## Success criteria

- `pip install -e .` on a documented CI image (gfortran + OpenMP + BLAS; no LAPACK or NetCDF Fortran needed for the extension).
- A driver script builds `PlasmaSurface` / `CoilSurface` / `Regcoil` and gets a `Solution` matching a legacy example within tolerance; `scan(...)` and `solve_for_target(...)` match `lambda_search_*`.
- Two `Regcoil` instances with different resolutions coexist and don't interfere (the kernel is stateless).
- `Regcoil` symmetry is controlled by the boolean `stellarator_symmetric`.
- Example regressions pass; pytest unit tests cover the Python object model and the Fortran kernels.
- `regcoil plot|compare|cut-coils` (names TBD) replace standalone scripts; Plotly coil visualization exists; no `.m` files remain.
- SVD scan, adjoint, WSO, and Laplace–Beltrami are gone; no `regcoil_variables`, handle API, or in-Fortran solve remain.
- Read the Docs publishes the user manual; LaTeX `manual/` retired; `publish_manual.yml` deleted.
- Dependency set matches the principles above (no SIMSOPT/DESC).

## Working agreements for agents

- Follow [PHASES.md](PHASES.md) exit criteria; stop when the phase is done.
- Prefer parity with examples before deleting the legacy Fortran driver or solve chain.
- Keep the Fortran boundary stateless and pure (ADR-020): explicit extents, `info` return, no `stop`, no module state, no LAPACK. Do all linear algebra in numpy/scipy (ADR-021).
- Do not add dependencies outside the allowed list without a new ADR.
- Use explicit typed options (for `Regcoil`, boolean `stellarator_symmetric`); parameters are set in Python code, not input files (ADR-019).
