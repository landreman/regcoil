# pyREGCOIL migration overview

Convert REGCOIL from a standalone Fortran executable into a **pip-installable Python package** that drives a slim, **instance-based** Fortran extension for the expensive linear algebra. Python owns I/O, orchestration, plotting/coil tools, and documentation; Fortran owns matrix assembly and LAPACK solves—with **no module-level globals** so multiple problems can coexist in memory.

Companion docs:

- [PHASES.md](PHASES.md) — ordered work packages and exit criteria
- [API.md](API.md) — intended Python / Fortran boundary
- [DECISIONS.md](DECISIONS.md) — append-only design log
- [INVENTORY.md](INVENTORY.md) — current-code map (roles, kill lists, touch points)

## Goals

1. **Python packaging** — `pip install` builds the Fortran extension and exposes a Python API / CLI.
2. **Python-driven execution** — remove the Fortran `PROGRAM regcoil` as the user-facing entry point.
3. **Keep Fortran where it matters** — matrix assembly, prepare/solve (BLAS/LAPACK), and related hot paths stay compiled.
4. **No Fortran globals** — replace `regcoil_variables` module state with derived-type / explicit arguments (and Python class state) so multiple instances can hold different values.
5. **Clean lambda search** — root-finding separated from physics; drive with SciPy.
6. **String-valued qualitative options** — integer “option codes” for qualitatively different modes become strings (less error-prone); see API.md.
7. **Dual input formats** — keep Fortran namelist input via **`f90nml`**, and also accept **JSON**.
8. **NetCDF only from Python** — no NetCDF compile/link dependency in Fortran.
9. **Fold tooling into the package** — `regcoilPlot`, `compareRegcoil`, and coil-cutting scripts become `pyregcoil` modules / CLI entry points.
10. **Plotly coil plot** — port `m20160811_01_plotCoilsFromRegcoil.m` to Python with Plotly; **remove all other MATLAB** sources.
11. **Remove SVD scan** — delete `regcoil_svd_scan` and `general_option == 3` pathways.
12. **Remove adjoint / winding-surface optimization** — delete sensitivity / adjoint Fortran and `windingSurfaceOptimization/`.
13. **Unit + regression tests** — pytest for Python and integration; unit tests for Fortran cores (framework TBD; see ADR-008). CI on GitHub Actions.
14. **Modern docs** — replace the LaTeX manual with a Read the Docs site (Sphinx or equivalent).
15. **Layout** — `src/`, `tests/`, `docs/`, `fortran/` (or similar) + `pyproject.toml`.

## Design principles

- **Minimal dependencies.**
  - **Python (allowed):** `numpy`, `scipy`, `matplotlib`, `f90nml`, `plotly`, and optionally `netCDF4` (see ADR-004). **Testing:** `pytest`.
  - **Fortran:** BLAS/LAPACK (and a Fortran compiler). No NetCDF Fortran libs after Phase 7.
  - **Do not** depend on **SIMSOPT**, **DESC**, or other large stellarator stacks.
- **Pytest** is the Python test runner (unit, integration, examples).
- Prefer **instance-friendly** APIs (classes / derived types) over process-global state.
- Prefer **string enums** for qualitative choices; keep integers for counts, indices, and resolutions.
- One phase per PR when practical; do not mix directory moves with physics refactors.
- Record non-obvious choices in [DECISIONS.md](DECISIONS.md).

## Non-goals (this migration)

- Rewriting Biot–Savart / inductive matrix physics in Python.
- Preserving binary-identical NetCDF files if metadata/order changes slightly; **example regression tolerances** are the contract.
- Supporting adjoint-based winding-surface design or SVD scan in the new package.
- Porting `coilMetricScripts/*.m`, `regcoil.m`, or `m20160811_02_compare2CoilsetsFromRegcoil.m` (delete with other MATLAB).
- Introducing SIMSOPT/DESC (or similar) as dependencies.

## Target architecture

```text
pip package: pyregcoil
│
├── Python
│   ├── CLI / RegcoilProblem / RegcoilResult
│   ├── input: f90nml namelist + JSON (string options)
│   ├── geometry / VMEC / nescin / bnorm loaders
│   ├── lambda scan + auto-λ (scipy.optimize)
│   ├── NetCDF I/O (library TBD: ADR-004)
│   ├── plot / compare / cut-coils tools (matplotlib + plotly)
│   └── docs → Read the Docs
│
└── Fortran extension (compiled; BLAS/LAPACK only; no NetCDF)
      ├── derived type (or equivalent) per problem instance — no globals
      ├── build matrices
      ├── prepare + solve (fixed λ)
      └── diagnostics for one λ
```

## Success criteria

- `pip install -e .` on a documented CI image (gfortran + BLAS/LAPACK; no Fortran NetCDF).
- Multiple `RegcoilProblem` instances can exist concurrently with different parameters/results.
- Qualitative options are strings in JSON and in the Python API; namelist accepts the same string values (legacy integer codes may be rejected or briefly translated—see ADR-009).
- Example regressions pass; pytest unit tests cover Python and Fortran-callable cores.
- `pyregcoil plot|compare|cut-coils` (names TBD) replace standalone scripts; Plotly coil visualization exists; no `.m` files remain.
- SVD scan, adjoint, and WSO are gone.
- Input works from both namelist (`f90nml`) and JSON.
- Read the Docs publishes the user manual; LaTeX `manual/` retired.
- Dependency set matches the principles above (no SIMSOPT/DESC).

## Working agreements for agents

- Follow [PHASES.md](PHASES.md) exit criteria; stop when the phase is done.
- Prefer parity with examples before deleting the Fortran driver.
- Do not add dependencies outside the allowed list without a new ADR.
- When touching options or input, prefer strings for qualitative modes and keep both namelist and JSON paths in sync.
