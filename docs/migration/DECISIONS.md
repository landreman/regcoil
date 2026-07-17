# Design decisions (ADR log)

Append-only. When you make a non-obvious choice, add a dated entry. Do not rewrite history; supersede with a new ADR.

Status values: `proposed` | `accepted` | `superseded` | `rejected`

---

## ADR-001: Package and import name

- **Date:** 2026-07-17
- **Status:** superseded
- **Context:** Repo is `pyREGCOIL`; Fortran project is REGCOIL; Python imports conventionally lowercase.
- **Options:** `pyregcoil` / `PyREGCOIL` / `regcoil`
- **Decision (tentative):** Distribution and import name `pyregcoil`.
- **Consequences:** Docs and CLI use `pyregcoil`; GitHub repo name can stay as-is.
- **Superseded by:** ADR-015.

---

## ADR-002: Build backend for the Fortran extension

- **Date:** 2026-07-17
- **Status:** accepted
- **Context:** Need pip-installable build with gfortran + BLAS/LAPACK; NetCDF must not be required after Phase 8.
- **Options:**
  1. `meson-python` + meson Fortran
  2. `scikit-build-core` + CMake
  3. Keep makefile and invoke it from a custom PEP 517 hook (least preferred long-term)
- **Decision:** Use **`meson-python` + Meson**. Binding style (f2py vs `iso_c_binding`) is deferred to Phase 4; this ADR only fixes the PEP 517 backend.
- **Consequences:** `pyproject.toml` uses `build-backend = "mesonpy"`; Phase 4 wires the Fortran extension in `meson.build`. Avoid custom makefile hooks for the supported install path. Legacy `makefile` remains until Phase 12 (ADR-006).

---

## ADR-003: SciPy root finder for auto-regularization

- **Date:** 2026-07-17
- **Status:** proposed
- **Context:** `regcoil_auto_regularization_solve.f90` embeds Brent; `regcoil_fzero.f` is SLATEC Dekker/Brent. Goal is separation of search vs physics.
- **Options:** `scipy.optimize.brentq`, `root_scalar(method="brentq"|"ridder"|"bisect")`, keep Fortran Brent behind a callback
- **Decision (tentative):** Python bracketing + `brentq`/`root_scalar` on `log(λ)` or `λ` (match current staging behavior closely enough for example tolerances).
- **Consequences:** Fortran exposes “solve at λ → residual metrics” only; example `lambda_search_*` become the acceptance tests.

---

## ADR-004: NetCDF library on the Python side

- **Date:** 2026-07-17
- **Status:** proposed (**needs explicit choice**)
- **Context:** Tests today use `scipy.io.netcdf_file`. Output is written via Fortran `ezcdf`. VMEC `wout` may be NetCDF. Allowed deps may include optional `netCDF4`. Minimal-dependency principle applies.
- **Options:**
  1. **`netCDF4`** — full NetCDF4/HDF5 features; common for VMEC `wout`; extra dependency (and system HDF5/NetCDF C libs in some installs).
  2. **`scipy.io.netcdf_file`** — already pulled in via SciPy; classic NetCDF3-style API; limited vs NetCDF4; historically used by `regcoilPlot` / tests / compare scripts.
  3. Hybrid: read/write classic outputs with SciPy; use `netCDF4` only if VMEC/NetCDF4 inputs require it.
- **Decision:** TBD — choose before Phase 8 implementation.
- **Tradeoffs to weigh:** dependency weight vs ability to read modern VMEC NetCDF; write compatibility with existing test variable names; CI image complexity.
- **Consequences:** Fortran NetCDF disappears in Phase 8 regardless of choice.

---

## ADR-005: Fate of `mini_libstell`

- **Date:** 2026-07-17
- **Status:** proposed
- **Context:** Bundled `mini_libstell/` supplies kinds, constants, VMEC wout read, ezcdf, etc. Full STELLOPT `libstell` is an alternate makefile path.
- **Options:**
  1. Vendor a trimmed mini_libstell inside `fortran/` (no ezcdf/NetCDF)
  2. Depend on external libstell (discouraged for pip; also conflicts with minimal deps)
  3. Replace wout/nescin reads entirely in Python and keep only kinds/constants in Fortran
- **Decision:** TBD; prefer (1) or (3) so `pip install` stays self-contained and avoids STELLOPT coupling.
- **Consequences:** Phase 8 must list which mini_libstell files remain.

---

## ADR-006: Dual executable during transition

- **Date:** 2026-07-17
- **Status:** accepted
- **Context:** Risk of silent physics drift while moving the driver to Python.
- **Decision:** Keep building the Fortran `regcoil` executable until Python-path example parity is proven in CI; remove in Phase 12.
- **Consequences:** Slightly larger build; clearer rollback; PRs should say which entry path CI uses.

---

## ADR-007: Adjoint / winding-surface optimization

- **Date:** 2026-07-17
- **Status:** accepted
- **Context:** Overhaul removes adjoint methods and winding-surface optimization.
- **Decision:** Delete in Phase 1; do not port to Python; related namelist keys error clearly if present.
- **Consequences:** Smaller extension; WSO scripts and adjoint manual chapter removed.

---

## ADR-008: How to unit-test Fortran cores

- **Date:** 2026-07-17
- **Status:** proposed (**needs explicit choice**)
- **Context:** Need unit tests for matrix/solve code, not only full-example regressions. Python testing standard is **pytest**. Minimal dependencies—avoid heavy Fortran test frameworks unless clearly worth it.
- **Options:**
  1. **pytest + Python extension** — call Fortran through the regcoil extension with small manufactured inputs; assert residuals/norms. No new Fortran test dep. Couples unit tests to bindings existing (Phase 4+).
  2. **Fortran framework** (e.g. [test-drive](https://github.com/fortran-lang/test-drive), pFUnit) — tests live next to `.f90`; can run before Python bindings exist; adds a Fortran test dependency/build path in CI.
  3. **Hybrid** — thin Fortran smoke tests early; prefer pytest-driven tests once the extension is stable.
- **Decision:** TBD before Phase 10 ramps up (a spike during Phase 4–5 is enough to choose).
- **Recommendation lean:** prefer (1) or (3) to keep deps minimal and one test runner (`pytest`) in CI.
- **Consequences:** CI workflow shape; where developers add numeric regression cases.

---

## ADR-009: String options and legacy integer codes

- **Date:** 2026-07-17
- **Status:** proposed
- **Context:** Integer option codes are error-prone; several fields are already strings (`regularization_term_option`, `target_option`). Namelist via `f90nml` and JSON should share one vocabulary.
- **Options:**
  1. Strings only in new inputs; reject integers.
  2. Accept legacy integers for one release with deprecation warnings, then reject.
  3. Accept both indefinitely (discouraged).
- **Decision (tentative):** Prefer (2) during migration, then (1). Exact string spellings documented in RTD and mirrored in JSON schema / examples.
- **Consequences:** Example `regcoil_in.*` files need updating; validate_input logic moves to Python.

---

## ADR-010: SVD scan removal

- **Date:** 2026-07-17
- **Status:** accepted
- **Context:** Overhaul explicitly removes SVD scan (`general_option == 3`, `regcoil_svd_scan.f90`).
- **Decision:** Delete in Phase 1; do not expose in Python API.
- **Consequences:** Slightly smaller code and docs; any old inputs requesting SVD scan should error clearly.

---

## ADR-011: Allowed dependencies

- **Date:** 2026-07-17
- **Status:** accepted
- **Context:** Keep the stack small and avoid coupling to other stellarator frameworks.
- **Decision:**
  - **Python runtime:** `numpy`, `scipy`, `matplotlib`, `f90nml`, `plotly`; optionally `netCDF4` if ADR-004 selects it.
  - **Python test/dev:** `pytest` (plus doc builder for RTD, e.g. Sphinx—docs-only).
  - **Fortran:** BLAS/LAPACK (+ compiler).
  - **Forbidden as package deps:** SIMSOPT, DESC, and similar large stacks.
- **Consequences:** Any new dependency needs a superseding ADR; tooling must be implemented in-tree with the allowed libs.

---

## ADR-012: Input formats (namelist + JSON)

- **Date:** 2026-07-17
- **Status:** accepted
- **Context:** Preserve familiar Fortran namelists while making structured input easier for scripts and tests.
- **Decision:** Support both: parse namelists with **`f90nml`** in Python; also accept JSON with the same logical keys. Fortran-side namelist read is not required on the supported path after Phase 6.
- **Consequences:** Single validation layer in Python; examples can migrate to JSON gradually; round-trip tests namelist ↔ dict ↔ JSON.

---

## ADR-013: MATLAB removal and Plotly port

- **Date:** 2026-07-17
- **Status:** accepted
- **Context:** Overhaul removes MATLAB; one useful visualization should move to Plotly.
- **Decision:** Port only `m20160811_01_plotCoilsFromRegcoil.m` to Python/Plotly; delete all other `.m` files (`regcoil.m`, `m20160811_02_*.m`, entire `coilMetricScripts/`). Fold existing Python scripts (`regcoilPlot`, `compareRegcoil`, coil cutters) into the package (matplotlib where they already use it).
- **Consequences:** Phase 9 delivers in-package CLIs; no MATLAB runtime needed.

---

## ADR-014: Documentation stack

- **Date:** 2026-07-17
- **Status:** accepted
- **Context:** LaTeX manual + `.github/workflows/publish_manual.yml` (build PDF, deploy to gh-pages) is dated relative to a Python package.
- **Decision:** Replace with **Read the Docs** (Sphinx preferred default; MkDocs acceptable if chosen in Phase 11). Retire `manual/` as canonical docs. **Delete** `publish_manual.yml` entirely—do not migrate that workflow to RTD.
- **Consequences:** New `docs/` structure; no gh-pages LaTeX publish job; doc deps are docs-only, not runtime.

---

## ADR-015: Package and import name is `regcoil`

- **Date:** 2026-07-17
- **Status:** accepted
- **Context:** ADR-001 tentatively chose `pyregcoil`; prefer the shorter name matching the Fortran code and paper.
- **Decision:** Distribution name, import package, and CLI entry-point prefix are all **`regcoil`**.
- **Consequences:** Layout is `src/regcoil/`; `import regcoil`; console scripts like `regcoil run|plot|…` (names TBD). Repo / migration-doc branding may still say pyREGCOIL.

---

## ADR-016: Phase 3 CI scope (build + smoke; examples later)

- **Date:** 2026-07-17
- **Status:** accepted
- **Context:** Phase 3 wants GHA on Ubuntu and macOS plus a pytest scaffold. Full `examples/` regressions take several minutes and still use the legacy executable/`runExamples.py` path. Branch-protection “required checks” were dropped from Phase 3 exit criteria.
- **Options:**
  1. Build + pytest smoke on both OS; defer example suite in CI
  2. Build + smoke on both OS; examples only on Ubuntu
  3. Full examples on both OS
- **Decision:** **(1)** for now. Workflow `.github/workflows/ci.yml` builds via `make MACHINE=github_ubuntu|github_macos`, installs the Python package, runs `pytest` (unit smoke). Example regressions remain `make test` locally until a follow-up enables them in GHA (toward pytest).
- **Consequences:** CI stays fast and unblocks packaging work; example parity in CI is an explicit later task (still Phase 3 follow-up or Phase 10 continuous work).

---

## ADR-017: Binding style is `iso_c_binding`

- **Date:** 2026-07-17
- **Status:** accepted
- **Context:** ADR-002 fixed meson-python but deferred f2py vs `iso_c_binding`. Phase 5 needs instance handles; Phase 4 still allows Fortran globals.
- **Options:**
  1. **f2py** — faster first wrap with globals; awkward for opaque instance state later
  2. **`iso_c_binding` + thin C Python extension** — more boilerplate now; natural for Phase 5 handles
- **Decision:** **(2)**. Module `regcoil._core` is a C extension calling `bind(C)` entry points in `fortran/regcoil_c_api.f90`.
- **Consequences:** No f2py in the build; NumPy array exchange can be added later without changing the binding strategy; Phase 5 can evolve the C API toward an opaque problem pointer.
