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
- **Status:** superseded
- **Context:** `regcoil_auto_regularization_solve.f90` embeds Brent; `regcoil_fzero.f` is SLATEC Dekker/Brent. Goal is separation of search vs physics.
- **Options:** `scipy.optimize.brentq`, `root_scalar(method="brentq"|"ridder"|"bisect")`, keep Fortran Brent behind a callback
- **Decision (tentative):** Python bracketing + `brentq`/`root_scalar` on `log(λ)` or `λ` (match current staging behavior closely enough for example tolerances).
- **Consequences:** Fortran exposes “solve at λ → residual metrics” only; example `lambda_search_*` become the acceptance tests.
- **Superseded by:** ADR-021. A generalized eigendecomposition makes the whole λ
  family closed-form, so `chi2`/`max_K` vs λ are smooth scalar functions; a root
  find is trivial and needs no Brent staging at all.

---

## ADR-004: NetCDF library on the Python side

- **Date:** 2026-07-17
- **Status:** resolved — split between reading and writing (see below); superseded on the output side by ADR-028.
- **Resolution:** **Reading** VMEC `wout` uses `scipy.io.netcdf_file` with a
  `netCDF4` fallback (the hybrid, option 3; see the Phase 6 note in ADR-023).
  **Writing** REGCOIL's own save/load output uses NetCDF-4 via `h5netcdf`
  (ADR-028), not `netCDF4` and not classic NetCDF-3.
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
- **Status:** superseded
- **Context:** Integer option codes are error-prone; several fields are already strings (`regularization_term_option`, `target_option`). Namelist via `f90nml` and JSON should share one vocabulary.
- **Options:**
  1. Strings only in new inputs; reject integers.
  2. Accept legacy integers for one release with deprecation warnings, then reject.
  3. Accept both indefinitely (discouraged).
- **Decision (tentative):** Prefer (2) during migration, then (1). Exact string spellings documented in RTD and mirrored in JSON schema / examples.
- **Consequences:** Example `regcoil_in.*` files need updating; validate_input logic moves to Python.
- **Superseded by:** ADR-019. The qualitative symmetry choice in `Regcoil`
  carries forward into the Python API as boolean `stellarator_symmetric`
  (`True` => stellarator-symmetric basis only, `False` => both sine and cosine);
  the namelist/JSON/integer-translation machinery is dropped with the input-file
  format.

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
- **Amended by:** ADR-019 (drop `f90nml` with the input-file format) and ADR-021
  (drop Fortran BLAS/LAPACK *for the solve*; numpy/scipy provide it). Fortran
  BLAS is not fully dropped: `regcoil_build_g_and_h` (ADR-020) keeps a DGEMM
  call by design, chosen over the `matmul` intrinsic for performance — Fortran
  LAPACK goes away (Phase 9), Fortran BLAS does not. The forbidden-deps decision
  (no SIMSOPT/DESC) stands.
- **Context:** Keep the stack small and avoid coupling to other stellarator frameworks.
- **Decision:**
  - **Python runtime:** `numpy`, `scipy`, `matplotlib`, `f90nml`, `plotly`; optionally `netCDF4` if ADR-004 selects it.
  - **Python test/dev:** `pytest` (plus doc builder for RTD, e.g. Sphinx—docs-only).
  - **Fortran:** BLAS (+ compiler + OpenMP); no LAPACK.
  - **Forbidden as package deps:** SIMSOPT, DESC, and similar large stacks.
- **Consequences:** Any new dependency needs a superseding ADR; tooling must be implemented in-tree with the allowed libs.

---

## ADR-012: Input formats (namelist + JSON)

- **Date:** 2026-07-17
- **Status:** superseded
- **Context:** Preserve familiar Fortran namelists while making structured input easier for scripts and tests.
- **Decision:** Support both: parse namelists with **`f90nml`** in Python; also accept JSON with the same logical keys. Fortran-side namelist read is not required on the supported path after Phase 6.
- **Consequences:** Single validation layer in Python; examples can migrate to JSON gradually; round-trip tests namelist ↔ dict ↔ JSON.
- **Superseded by:** ADR-019. There is no input-file format; parameters are set in
  code via the object model. `f90nml` is dropped from the dependency list.

---

## ADR-013: MATLAB removal and Plotly port

- **Date:** 2026-07-17
- **Status:** accepted
- **Context:** Overhaul removes MATLAB; one useful visualization should move to Plotly.
- **Decision:** Port only `m20160811_01_plotCoilsFromRegcoil.m` to Python/Plotly; delete all other `.m` files (`regcoil.m`, `m20160811_02_*.m`, entire `coilMetricScripts/`). Fold existing Python scripts (`regcoilPlot`, `compareRegcoil`, coil cutters) into the package (matplotlib where they already use it).
- **Consequences:** Phase 10 delivers in-package plotting/compare/cut and a single `regcoil` CLI (architecture fixed in ADR-029); no MATLAB runtime needed.

---

## ADR-014: Documentation stack

- **Date:** 2026-07-17
- **Status:** superseded
- **Context:** LaTeX manual + `.github/workflows/publish_manual.yml` (build PDF, deploy to gh-pages) is dated relative to a Python package.
- **Decision:** Replace with **Read the Docs** (Sphinx preferred default; MkDocs acceptable if chosen in Phase 11). Retire `manual/` as canonical docs. **Delete** `publish_manual.yml` entirely—do not migrate that workflow to RTD.
- **Consequences:** New `docs/` structure; no gh-pages LaTeX publish job; doc deps are docs-only, not runtime.
- **Superseded by:** ADR-030 (keeps Read the Docs + retire `manual/` + delete `publish_manual.yml`; fixes the remaining open choices — Sphinx over MkDocs, theme, executed-notebook mechanism, and where code runs).

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

---

## ADR-018: Instance state layout and Fortran calling pattern (Phase 5)

- **Date:** 2026-07-17
- **Status:** superseded
- **Superseded by:** ADR-020. Phase 5 wrapped the in-Fortran solve behind an
  opaque `type(regcoil_t)` handle; the new architecture moves the solve to
  numpy/scipy and narrows Fortran to stateless pure functions, so no problem
  state (and no derived type) crosses the boundary. The handle API and
  `regcoil_variables.f90` are removed in Phases 7–9 / 13.
- **Context:** Phase 5 removes `regcoil_variables` module globals so multiple problems can coexist. Nested types were preferred over a single flat derived type.
- **Decision:**
  1. **`type(regcoil_t)`** holds:
     - `plasma` — `regcoil_plasma_surface_t`
     - `coil` — `regcoil_coil_surface_t`
     - `input` — `regcoil_solver_input_t` (namelist/solver parameters, including `lambda`)
     - `output` — `regcoil_solver_output_t` (diagnostics / potentials)
     - `work` — `regcoil_work_t` (matrices, basis, LAPACK scratch; needed companion to the four surface/solver types)
  2. **Calling convention:** new and updated Fortran routines take `type(regcoil_t), intent(inout) :: prob` as the **first** argument. Prefer `associate` for readability; do **not** `allocate`/`deallocate` a selector while it is associated—allocate first, then associate, or use `prob%…` for allocate targets.
  3. **fzero scratch** in offset-surface / plasma-init is **local** (nested residual via host association), not module variables.
  4. **mini_libstell** module state may remain until a later phase.
  5. **Legacy `program regcoil`** holds one local `regcoil_t` and passes it through (no compatibility globals).
  6. **Python:** opaque handle via `regcoil_c_create` / `destroy`; thin `RegcoilProblem` wraps setup/solve. Namelist still filled in Fortran for now (Phase 6 moves I/O).
- **Consequences:** C API is handle-based; concurrent-instance pytest is required; string option PARAMs stay as module constants (immutable).

---

## ADR-019: Scriptable object-model API (no input-file format)

- **Date:** 2026-07-18
- **Status:** accepted
- **Supersedes:** ADR-012 (namelist + JSON input); ADR-009 (string-enum decision
  carried forward, integer/namelist machinery dropped).
- **Context:** The overhaul (design session 2026-07-18) drops the input-file
  format entirely: users drive the calculation with Python scripts and set
  parameters in code. The two natural objects are the plasma surface and the coil
  surface; existing initialization methods must survive and stay extensible.
- **Decision:** User-facing objects are `PlasmaSurface`, `CoilSurface`, `Regcoil`,
  and a `Solution` result, over a `Surface` (ABC) / `FourierSurface` base. Legacy
  `geometry_option_*` codes become **alternate constructors** (`from_wout`,
  `from_nescin`, `from_focus`,`from_ascii_table`, `circular_torus`,
  `from_uniform_offset`). The `Surface` ABC (`_evaluate` contract) is the
  extensibility hook. `Regcoil` accepts a `Surface` (documented arrays it reads
  let an outer optimization loop supply its own geometry). Qualitative options are
  **string enums** in the API. No namelist/JSON reader; `f90nml` leaves the
  dependency list.
- **Consequences:** See [API.md](API.md). `regcoil_read_input`/`write_input`/
  `validate_input` and the Python `RegcoilProblem`/`RegcoilResult` sketch are
  retired; docs (Phase 12) document the object model, not input parameters.

---

## ADR-020: Stateless Fortran kernel boundary

- **Date:** 2026-07-18
- **Status:** accepted
- **Supersedes:** ADR-018 (opaque instance handle / derived types crossing the
  boundary / in-Fortran solve).
- **Context:** Profiling the design (attached session) shows only two genuinely
  expensive hot paths that belong in Fortran: the O(N²) pairwise inductance/field
  sum, and the uniform-offset surface (serial root solves + splines + DFT).
  Everything else is a gemm or a small solve that numpy/scipy already do with
  threaded BLAS/LAPACK.
- **Decision:** Delete `regcoil_variables.f90`. Fortran is **pure** subroutines
  with `intent(in)`/`intent(out)` arrays and **explicit** extents: no module
  globals, no derived types crossing the boundary, no persistent handle, **no
  `stop`** (every entry returns `integer :: info`), and **no LAPACK**. Entry
  points: `regcoil_build_g_and_h` (fused `inductance @ basis`, OpenMP,
  threadsafe), `regcoil_build_inductance` (debug/regression, full matrix), and
  `regcoil_uniform_offset_surface` (returns Fourier coefficients). Optional
  `regcoil_evaluate_surface` only if profiling demands. `iso_c_binding` +
  `regcoil._core` (ADR-017) and `meson-python` (ADR-002) are unchanged — only the
  number/shape of entry points changes.
- **Consequences:** The extension links only the compiler runtime, OpenMP, and
  BLAS after Phase 9 (BLAS stays for `regcoil_build_g_and_h`'s DGEMM, chosen
  over the `matmul` intrinsic for performance; LAPACK and NetCDF are the ones
  that go away — Phase 9 retitled accordingly). Kernels are unit-testable in
  isolation from Python. The handle API and in-Fortran solve/prepare/diagnostics
  are removed across Phases 7–9 / 13.

---

## ADR-021: Linear algebra and λ family in numpy/scipy

- **Date:** 2026-07-18
- **Status:** accepted
- **Supersedes:** ADR-003 (Fortran → SciPy Brent auto-regularization).
- **Context:** numpy's `@`/`dot` and `scipy.linalg` call the same threaded
  BLAS/LAPACK a Fortran layer would; reimplementing them in Fortran buys nothing
  and costs an interface. `matrix_B` and `matrix_K` are symmetric with `matrix_K`
  positive definite.
- **Decision:** Assemble `matrix_B`, `matrix_K`, and the RHS vectors in numpy, and
  solve in scipy. Compute **one** generalized eigendecomposition
  `w, V = scipy.linalg.eigh(matrix_B, matrix_K)` in `Regcoil.__init__`; then
  `x(λ) = V @ ((Vᵀ·b_B + λ·Vᵀ·b_K) / (w + λ))` gives every λ in O(nbf²), and
  `f_B(λ)`, `f_K(λ)`, `max_K(λ)` are closed-form for scans and target
  searches. Retire `regcoil_auto_regularization_solve.f90`,
  `regcoil_lambda_scan.f90`, `regcoil_compute_lambda.f90`, and (for λ) `regcoil_fzero.f`.
- **Consequences:** No Fortran BLAS/LAPACK link *for the solve* (the unrelated
  DGEMM call inside `regcoil_build_g_and_h`, ADR-020, stays); the λ scan and
  target search are ~15 lines of Python; the object is immutable after
  construction (the expensive step is exactly the geometry-dependent one, so
  there is nothing to cache).

---

## ADR-022: Remove Laplace–Beltrami regularization

- **Date:** 2026-07-18
- **Status:** accepted
- **Context:** Open decision from the design session; the adjoint method is already
  being removed (ADR-007). Laplace–Beltrami is the only consumer of surface second
  derivatives and the metric-coefficient block in `build_matrices`.
- **Decision:** Drop Laplace–Beltrami regularization. Surface `_evaluate` returns
  first derivatives only (no `nderiv`), and the second-derivative / metric block
  disappears from matrix assembly.
- **Consequences:** Simpler `Surface` contract and matrix build; the
  `regularization_term_option` choices reduce accordingly.

---

## ADR-023: Phase 6 scope decisions (surface object model)

- **Date:** 2026-07-19
- **Status:** accepted
- **Context:** Implementing `Surface`/`FourierSurface`/`PlasmaSurface`/`CoilSurface`
  (Phase 6) surfaced three items not fully pinned down by ADR-019/API.md.
- **Decisions:**
  1. **`from_efit` dropped, not stubbed-and-deferred.** Legacy Fortran explicitly
     stops with "geometry_option_plasma=5 (EFIT) is no longer supported"
     (`regcoil_validate_input.f90`); there is no reference implementation or test
     fixture to build or validate a new EFIT g-file reader against. `from_efit` is
     omitted from `PlasmaSurface` and from the API.md constructor table until a
     concrete need (and a g-file + reference output to validate against) exists.
  2. **`from_wout(straight_field_line=True)` raises `NotImplementedError`.** The
     legacy coordinate transform (geometry_option_plasma=4) root-solves VMEC's
     theta against a fixed ±0.3 bracket per grid point; this is not robust even in
     the reference Fortran (observed "no sign change in residual" errors on a
     real VMEC file at moderate transform resolution). Porting it without a
     validated reference risked shipping a silently-wrong surface. Deferred past
     Phase 6; `mesh="full"`/`mesh="half"` are implemented and legacy-matched.
  3. **VMEC `wout` reading uses `scipy.io.netcdf_file` first, `netCDF4` as a
     fallback** for HDF5-based wout files (`src/regcoil/_io.py:read_vmec_wout`).
     This is the "hybrid" option from ADR-004, chosen because scipy is already a
     required dependency and reads both example `wout` files; `netCDF4` stays
     optional. ADR-004 remains open for Phase 9 output writing.
- **Consequences:** `CoilSurface.from_uniform_offset` (Phase 7, needs the Fortran
  offset kernel) and `from_efit` both raise/omit rather than half-implement;
  `PlasmaSurface.from_wout`, `from_ascii_table`, `from_focus`,
  `CoilSurface.from_nescin`, and `FourierSurface.circular_torus` are implemented
  and checked against the legacy Fortran (`regcoil_init_plasma`/
  `regcoil_init_coil_surface`, compiled standalone for comparison) in
  `tests/unit/`.

---

## ADR-024: Phase 8 scope decisions (`Regcoil`/`Solution` assembly and solve)

- **Date:** 2026-07-19
- **Status:** accepted
- **Context:** Implementing `Regcoil`/`Solution` (Phase 8) surfaced a few items
  not fully pinned down by ADR-021/API.md.
- **Decisions:**
  1. **`regularization_term_option` is f_K only.** `matrix_K`/`RHS_K` are
     built from `f_x, f_y, f_z` (f_K's definition) only; the legacy `K_xy`
     and `K_zeta` regularization variants (and Laplace–Beltrami, already
     dropped by ADR-022) are not ported -- no example or ADR calls for them,
     and API.md's `matrix_K = Σ fᵢᵀ(fᵢ/N)` already specifies exactly this.
  2. **`solve_for_target` raises `ValueError` for an unreachable target**
     (the achievable range is read off the `lam=0`/`lam=inf` closed-form
     endpoints) instead of the legacy `exit_code` sentinel
     (`-2`/`-3`/other). There is no `RegcoilProblem`/`RegcoilResult` handle to
     stash a sentinel on post-ADR-019, and a Python API should signal
     "can't do that" via an exception, not a magic return code.
  3. **`solve_for_target` bisects in `log(lambda)`** via `scipy.optimize.brentq`
     rather than porting `regcoil_auto_regularization_solve.f90`'s staged
     bracket-then-Brent search. The lambda *sequence* visited necessarily
     differs from the legacy solver (already anticipated in PHASES.md Phase
     8); only the converged endpoint is a contract.
  4. **`CoilSurface.from_uniform_offset` stands in for legacy
     `geometry_option_coil=4`** (the constant-arclength theta
     reparametrization) in the `regcoilPaper_figure10d_but_geometry_option_coil_4_loRes`
     regression test. The Phase 7 kernel
     (`regcoil_uniform_offset_surface`) only ports plain `geometry_option_coil=2`
     (ADR-020); there is no Python constructor for the constant-arclength
     iteration, and no ADR calls for adding one. This is legitimate because the
     legacy example's own point is that the physics is independent of the
     coil-surface parametrization -- its golden chi2/max_K tolerances (3-6%)
     already budget for exactly this kind of small parametrization difference.
  5. **`Solution.single_valued_current_potential_mn` is an alias for
     `Solution.solution`** (a `@property`, not a stored duplicate) -- in the
     legacy Fortran the two are the same array (`single_valued_current_potential_mn(:,ilambda)
     = solution`); API.md lists both names only because it mirrors the legacy
     NetCDF variable name for Phase 9's `Solution.save`.
- **Consequences:** Regression tests in `tests/regression/` build the problem
  directly via `PlasmaSurface`/`CoilSurface`/`Regcoil` (no legacy executable,
  no NetCDF) and compare against golden values read (or, for the large
  full-resolution mode-coefficient arrays, extracted programmatically) from
  `examples/*/tests.py`; `ntheta_plasma=128` cases are `@pytest.mark.slow`
  and skipped in CI (`pytest -m "not slow"`) per PHASES.md Phase 8.

---

## ADR-025: `Surface.standard_toroidal_angle`; second `from_uniform_offset` construction

- **Date:** 2026-07-19
- **Status:** accepted
- **Context:** `CoilSurface.from_uniform_offset` (Phase 7) only implemented the
  legacy `geometry_option_coil=2` construction: for each coil `(theta, zeta)`
  grid point, root-solve (`regcoil_uniform_offset_surface`) for the plasma
  point whose normal-offset lands at Cartesian toroidal angle `atan2(y,x) ==
  zeta`. Future cross-section plotting (Phase 10) plots a surface's `r` at
  fixed `zeta`, which is only a plane of constant physical toroidal angle if
  `zeta` *is* the standard toroidal angle in that sense. A second, simpler
  offset construction is wanted: move each plasma grid point `separation`
  along its local unit normal and Fourier-transform the result labeled by the
  *same* `(theta, zeta)`, with no root solve. Because the plasma normal
  generally has a toroidal component, this mislabels the moved point's actual
  `atan2(y,x)` as `zeta` -- a real geometric difference from the root-solve
  construction, not just an implementation detail, and one that would
  silently corrupt any future cross-section plot that assumes constant-`zeta`
  == constant-physical-angle.
- **Decision:**
  1. Add `Surface.standard_toroidal_angle: bool` (alongside `nfp`,
     `stellarator_symmetric`, ...): `True` iff a constant-`zeta` slice of `r`
     is a constant-physical-toroidal-angle cross section. All existing
     constructors (`circular_torus`, `from_wout`, `from_ascii_table`,
     `from_focus`, `from_nescin`, and the root-solve `from_uniform_offset`)
     set/default it `True`.
  2. `CoilSurface.from_uniform_offset` gains a `standard_toroidal_angle: bool`
     argument that both selects the construction and is stored as the
     result's attribute of the same name, so the argument name *is* the
     contract. **Default is now `False`** (the new normal-offset method,
     pure Python/numpy, no Fortran call) rather than the legacy root-solve;
     `standard_toroidal_angle=True` keeps the old behavior (and is required
     to reproduce the legacy-Fortran-golden regression tests, which pin it
     explicitly). `ntheta_transform`/`nzeta_transform` apply to both
     constructions; `tol` (root-solve tolerance) only affects `True`.
  3. The normal-offset method samples one field period, not the full torus:
     the plasma's `xn` is already an nfp multiple, so a `2*pi/nfp` rotation
     maps the plasma surface (and hence its normal field and any
     fixed-separation offset of it) to itself, and one period is enough to
     reconstruct the Fourier coefficients -- the same shortcut the root-solve
     kernel already relies on.
- **Consequences:** For an axisymmetric plasma the two constructions agree
  exactly (the normal has no toroidal component), giving a cheap
  cross-check (`tests/unit/test_coil_surface.py::test_from_uniform_offset_circular_torus_is_exact`,
  parametrized over both). For a non-axisymmetric plasma they differ by an
  amount that shrinks with `separation` and with how non-circular the
  cross-section is; both are physically legitimate offset surfaces (as
  ADR-024 item 4 already established for `=2` vs `=4`), just parametrized
  differently. All pre-existing unit/regression call sites that depend on
  matching the legacy Fortran (`tests/unit/test_kernels.py`,
  `tests/regression/lambda_search_*`,
  `tests/regression/regcoilPaper_figure10d_but_geometry_option_coil_*`) were
  updated to pass `standard_toroidal_angle=True` explicitly. Phase 10
  cross-section plotting must check this attribute (and, if plotting a
  `standard_toroidal_angle=False` surface, either say so or resample/re-root-
  solve rather than assume constant-`zeta` is a physical cross section).
- **Refined by:** ADR-026. ADR-025's first implementation of the
  `standard_toroidal_angle=False` construction was geometrically wrong (it
  DFT-ed the moved points' `(R, Z)` into a plain `FourierSurface`, which
  reconstructs at `phi = zeta` and so relocated every point toroidally); the
  attribute/argument/field-period decisions above stand, but the surface
  representation is fixed in ADR-026.

---

## ADR-026: `nu` toroidal-angle-shift Fourier modes for non-standard-toroidal-angle surfaces

- **Date:** 2026-07-20
- **Status:** accepted
- **Context:** ADR-025's `standard_toroidal_angle=False` construction moved each
  plasma point `separation` along its normal, took `major_R = sqrt(x²+y²)` and
  `Z = z` of the moved point, DFT-ed those as functions of `(theta, zeta)`, and
  built a plain `FourierSurface`. But `FourierSurface._evaluate` hard-codes the
  physical toroidal angle `phi = zeta` (`X = R*cos(zeta)`, `Y = R*sin(zeta)`).
  The moved point's true angle is `phi = atan2(y,x) = zeta + delta` with
  `delta != 0` whenever the plasma normal has a toroidal component, so the
  reconstructed surface was the offset surface with every point *rotated back*
  toroidally to `phi = zeta` -- a genuinely different (wrong) surface. The
  original tests only exercised axisymmetric plasmas, where `delta == 0` hides
  the bug.
- **Key insight:** representing a "moved along the normal" surface in Fourier
  space is only possible if the toroidal-angle deviation is carried explicitly.
  DFT-ing `(R, Z)` *without* it is exactly what broke. So either keep the
  surface as an exact geometric object (no Fourier rep) or add angle-shift
  modes. We chose the Fourier route (Option 1 from the design discussion),
  because a coil winding surface benefits from being a smooth, resolution-
  controlled, `filter_modes`-able Fourier surface (feature parity with the
  `standard_toroidal_angle=True`/legacy path), and it stays within the
  first-derivative-only surface contract (ADR-022) with no live plasma
  reference.
- **Decision:**
  1. Generalize `FourierSurface` with a **toroidal-angle-shift** mode set `nu`
     (`numns` sine part, `numnc` cosine part; same `(xm, xn)` mode numbering
     and `m*theta - n*zeta` angle as `R`/`Z`). The Cartesian point is placed at
     `phi = zeta + nu`: `x = R*cos(zeta+nu)`, `y = R*sin(zeta+nu)`, `z = Z`.
     Derivatives stay analytic (`dphi/dtheta = dnu/dtheta`,
     `dphi/dzeta = 1 + dnu/dzeta`) -- no second derivatives. Default `nu = 0`
     reduces exactly to the previous behavior (a fast path keeps the common
     `phi = zeta` case unchanged), so every existing standard-angle surface is
     byte-for-byte unaffected.
  2. Under stellarator symmetry `nu` has **sine parity** (`numnc = 0`), because
     the stellarator-symmetry map sends `phi -> -phi` and `phi = zeta + nu`, so
     `nu` must be odd -- the same parity as `Z`'s `zmns`. `numnc != 0` therefore
     also breaks the `stellarator_symmetric` auto-inference.
  3. `standard_toroidal_angle` stays an **explicit constructor argument**, not
     derived from whether `nu` is (numerically) zero. Deriving it is fragile:
     an axisymmetric offset yields `nu` at rounding-noise level, and a
     construction-intent flag is what downstream code (Phase 10 plotting) wants
     anyway. `from_uniform_offset` sets the attribute from its argument and
     supplies `nu` modes only on the `False` path; the two are always
     consistent because that constructor controls both.
  4. `from_uniform_offset`'s `False` path computes `nu = atan2(y,x) - zeta`
     (wrapped into `(-pi, pi]`; the true shift is `~separation/R << pi`, so the
     branch is unambiguous) on the field-period grid and DFTs it alongside
     `R`/`Z`. `filter_modes` zeroes `numns`/`numnc` with `R`/`Z` so smoothing
     stays consistent.
- **Consequences:** With `mpol`/`ntor` at the transform grid's Nyquist limit
  the DFT->IDFT round-trips, so the reconstructed coil surface passes through
  the moved-along-normal points to machine precision -- verified on a
  non-axisymmetric plasma in
  `tests/unit/test_coil_surface.py::test_from_uniform_offset_reproduces_moved_points_nonaxisymmetric`
  (this test fails on the pre-ADR-026 code). Analytic `nu` derivatives are
  checked against finite differences
  (`tests/unit/test_surface.py::test_nu_angle_shift_places_point_at_zeta_plus_nu`),
  so `Regcoil`'s coil normal and basis functions are correct on these
  surfaces. NetCDF output (Phase 9) and any nescin export must now carry the
  `nu` modes to round-trip a `standard_toroidal_angle=False` surface; a plain
  nescin file (R/Z only) cannot represent one without loss. The
  `standard_toroidal_angle=True` (Fortran root-solve) path is unchanged and
  produces `nu = 0`.

---

## ADR-027: Rename `Solution.chi2_B`/`chi2_K` to `f_B`/`f_K`

- **Date:** 2026-07-20
- **Status:** accepted
- **Context:** The legacy Fortran names `chi2_B`/`chi2_K` are misleading in the
  Python API: they are not chi-squared statistics, just the surface-integrated
  squared-field/current-density objectives being traded off by `lambda`.
- **Decision:** Rename the `Solution` dataclass fields (and the local variables
  that build them in `Regcoil._build_solution`) to `f_B`/`f_K`. Applies
  throughout the Python package -- `src/regcoil/`, `tests/unit/`,
  `tests/regression/` (including golden-value module names `CHI2_B`/`CHI2_K`
  -> `F_B`/`F_K`) -- and this migration documentation. The legacy Fortran
  source, the `examples/*/tests.py` netCDF-variable checks (which read the
  Fortran output variables literally named `chi2_B`/`chi2_K`), and the manual
  are untouched, since they describe the Fortran code being replaced.
- **Consequences:** Any external code calling `solve_for_target("chi2_B", ...)`
  or reading `Solution.chi2_B`/`chi2_K` must switch to `"f_B"`/`f_B`/`f_K`.

---

## ADR-028: Save/load format and object serialization (Phase 9)

- **Date:** 2026-07-20
- **Status:** accepted
- **Amended:** 2026-07-20 (see ADR-029) — the `/problem` group additionally
  stores **`Bnormal_from_net_coil_currents`** so the full Bnormal plot set is
  reproducible from a saved file with no Fortran kernel call. See decision 2
  below and the schema in [API.md](API.md#file-layout).
- **Supersedes:** ADR-004 (NetCDF library on the Python side) — the output-side
  question ADR-004 left open for Phase 9 is resolved here. VMEC `wout` *reading*
  stays as ADR-004 decided (`scipy.io.netcdf_file`, `netCDF4` fallback; see the
  Phase 6 note in ADR-023); this ADR governs REGCOIL's own saved output.
- **Context:** Phase 9's intent is object serialization — saving and loading
  `PlasmaSurface`, `CoilSurface`, `Regcoil`, and `Solution` objects (and the
  `Solution` list from a lambda scan), typically all together but optionally one
  at a time. There is **no requirement to preserve the legacy Fortran
  `regcoil_out.*.nc` format**; this is a fresh Python-defined schema. The driving
  use case is that everything a user might want to plot from a saved run — 3D
  surface shapes, surface cross-sections, the Pareto front between `f_B`/`f_K` or
  `max_Bnormal`/`max_K`, the current potential and current density on the
  `(theta, zeta)` grid, `Bnormal` on the plasma — must be available quickly
  (order 1 s at 128×128) **without calling any Fortran kernel or expensive BLAS,
  and ideally without reconstructing any Python object at all.**
- **Decisions:**
  1. **Format: NetCDF-4 written via `h5netcdf`.** HDF5 underneath (no NetCDF C
     library needed in wheels/CI), but the NetCDF-4 data model gives named,
     shared dimensions (`mnmax_plasma`, `mnmax_coil`, `nbf`, `nlambda`,
     `ntheta_coil`, …) so the file is self-documenting and directly
     `xarray.open_dataset`-able — a user can scatter `f_B` vs `f_K` for a Pareto
     front with no `regcoil` import. `h5netcdf` becomes a **required** dependency
     for I/O (core solve stays independent of it). Rejected: raw `h5py` (throws
     away free named dimensions for no benefit here); `netCDF4` (equivalent data
     model but drags in the NetCDF/HDF5 C libraries); `scipy.io.netcdf_file`
     (NetCDF-3: no groups, no compression — can't express the object graph);
     JSON (bloats/precision-risks the large float arrays). The ecosystem-tradition
     argument for NetCDF (VMEC `wout`, legacy `regcoil_out`) is explicitly **not**
     weighted — HDF5 is equally common in the fusion ecosystem.
  2. **Store the source of truth, plus exactly those derived quantities whose
     recomputation needs the Fortran kernel or the big operators.** The dividing
     line: recompute anything derivable from just the surfaces + the mode
     amplitudes; store anything that would otherwise require `g`, `h`,
     `matrix_B`/`matrix_K`, or the eigendecomposition `V`. **The big operators
     are never saved** (by default or otherwise in Phase 9). Concretely:
     - **Surfaces** (`PlasmaSurface`/`CoilSurface`): store only source of truth —
       `xm, xn, rmnc, rmns, zmnc, zmns, numns, numnc` (the `nu` modes per
       ADR-026), `nfp, ntheta, nzeta, stellarator_symmetric,
       standard_toroidal_angle`. `r`, `normal`, `area`, `volume` are one gemm
       away, recomputed on load. `PlasmaSurface` additionally stores
       `net_poloidal_current`, `curpol`, and **`Bnormal_from_plasma_current`**
       (the last is unrecoverable — it is set after construction from a
       bnorm/FOCUS file or a user array).
     - **`Regcoil`**: store scalar params (`mpol_potential, ntor_potential,
       net_poloidal_current, net_toroidal_current, stellarator_symmetric`), a
       reference to its plasma and coil, **and `Bnormal_from_net_coil_currents`**
       (the per-problem, λ-independent grid `h`-reshaped/`|N|`; it needs the
       kernel output `h` to recompute, so it is stored rather than rebuilt --
       added by the Phase-10 amendment below so the
       `bnormal(component="net_coil")` plot panel loads without a kernel call).
       `basis_functions`, `f_all`, `d_xyz` are cheap numpy and recomputed on
       load; `g, h, matrix_B/K, RHS_*, w, V` are **not** stored. **Exception (ADR-029):** the plasma-grid array
       `Bnormal_from_net_coil_currents(ntheta_plasma, nzeta_plasma)` *is* stored
       in `/problem`. It is λ-independent, ~0.1 MB, and its recompute is the
       Fortran `h` term (`regcoil_build_g_and_h`); storing it lets the Bnormal
       plot panel "from net coil currents" render from a saved file with no
       kernel call, matching the standard set every other plot target already
       satisfies. It is the only operator-derived quantity stored on `/problem`.
     - **`Solution`**: store `lam`, `solution` (the single-valued current-potential
       Fourier amplitudes — source of truth), the scalar diagnostics `f_B, f_K,
       max_K, rms_K, max_Bnormal`, **and the result-sized grids `Bnormal_total`,
       `current_potential`, and `current_density`.**
  3. **Store the result-sized grids rather than recompute them.** `Bnormal_total`
     *must* be stored (its recompute is `g @ solution`, needing the expensive
     kernel). `current_potential`/`current_density` *could* be recomputed from
     the amplitudes, but doing so materializes the `f_all` intermediate
     (`(nbf, 3, ncoil_grid)` ≈ 0.5 GB and ~1 s at 128×128 with `nbf ≈ 1250`),
     whereas the *results* are tiny (~0.4 MB per λ for all three grids; ~26 MB
     for a 40-λ scan). Storing the small results gives zero-compute,
     zero-object-reconstruction plotting and avoids the memory spike, at the cost
     of harmless redundancy with the stored amplitudes. This is the mechanism by
     which the "quick plotting, nothing expensive" requirement is met.
  4. **A lambda scan is flattened on disk along a `lambda` dimension**, not stored
     as per-solution groups. A scan is homogeneous by construction (one shared
     `Regcoil` → one plasma, one coil, one `nbf`, one grid), and the dominant
     analysis reads (`f_B(:)`, `max_K(:)`, … for a Pareto front) want columnar
     arrays. One `/solutions` group holds `lam(nlambda)`, `solution(nlambda,
     nbf)`, `f_B(nlambda)`, …, `Bnormal_total(nlambda, ntheta_plasma,
     nzeta_plasma)`, etc. A single `solve()` result is the `nlambda = 1` case.
     This flattening is chosen on its own merits, **not** because the legacy
     Fortran indexed by `ilambda`.
  5. **In memory, a scan is a `SolutionScan`** — a `Sequence[Solution]` that
     iterates/indexes as `Solution` objects (uniform with `solve()`'s return and
     preserving each solution's lazy `current_potential()`/`current_density()`)
     **and** exposes vectorized `.lam`, `.f_B`, `.f_K`, `.max_K`, `.rms_K`,
     `.max_Bnormal` arrays for direct Pareto/scan plotting. `Regcoil.scan()` and
     `regcoil.load()` both return it (a change from `scan()`'s current plain
     `list`).
  6. **Object graph and API.** File layout is groups by object with a `_class`
     attribute per group and a `format_version` at the root (for dispatch and
     future migration): `/plasma`, `/coil`, `/problem`, `/solutions`. The API is
     a polymorphic module-level pair plus thin instance methods:
     - `regcoil.save(path, *, plasma=None, coil=None, problem=None,
       solutions=None)` writes only the groups it is given, deduping the shared
       problem/surfaces (never N copies of the problem for an N-λ scan).
       `solutions` transitively carries its `problem` → `coil`, `plasma`, so
       `regcoil.save(path, solutions=scan)` writes the whole bundle.
     - `obj.save(path)` delegates to `regcoil.save` with just that object.
     - `regcoil.load(path)` returns a container exposing `.plasma`, `.coil`,
       `.problem`, `.solutions` (any may be `None`); `Class.load(path)` returns
       the single object, erroring if that group is absent.
  7. **Loaded objects and the cheap/expensive assembly split.** To make a loaded
     `Regcoil` power its Solutions' `current_potential()`/`current_density()`
     without the operators, `Regcoil.__init__` is split into a **cheap** assembly
     (surfaces, potential modes, `basis_functions`, `f_all`, `d_xyz` — no Fortran,
     no big BLAS) and an **expensive** assembly (`build_g_and_h`, `matrix_B/K`,
     RHS, `eigh`). `load()` runs only the cheap part; a loaded `Regcoil` gains
     operators lazily (re-running the expensive assembly, with a log message) only
     if the user asks for a *new* λ via `solve()`/`scan()`. Reconstructing already
     saved Solutions for plotting never triggers this.
- **Consequences:**
  - `h5netcdf` is added to the runtime dependencies in `pyproject.toml`
    (updating the ADR-011 dependency list); the core solve remains importable and
    runnable without it, but `save`/`load` require it.
  - Saved files are readable by generic tooling (`xarray`, `h5py`,
    `ncdump`) with no `regcoil` import, satisfying the "plot from saved output"
    goal directly; and by `regcoil.load` for full object reconstruction.
  - Redundancy between the stored result grids and the stored amplitudes is
    accepted; the file is a data product, not a minimal object pickle, and a
    hand-edited file could make the two inconsistent.
  - A `standard_toroidal_angle=False` coil surface round-trips only because the
    `nu` modes are in the stored source of truth (ADR-026); there is no lossy
    R/Z-only path.
  - PHASES.md Phase 9 is rewritten around "save/load object serialization" (the
    old "strip Fortran NetCDF/LAPACK" bullets stay as the parallel cleanup they
    always were, now clearly secondary to the serialization work). Phase 11's
    "NetCDF round-trip" test becomes a round-trip test against this schema.
- **Amendments:**
  - *(Phase 10, ADR-029):* `/problem` also stores `Bnormal_from_net_coil_currents`
    (the λ-independent, `h`-derived grid) so the `bnormal(component="net_coil")`
    plot panel loads without a kernel call. Folded into decision 2's `Regcoil`
    storage list above.

---

## ADR-029: Plotting and visualization architecture (Phase 10)

- **Date:** 2026-07-20
- **Status:** accepted
- **Amends:** ADR-013 (the "fold the scripts in" intent stands; this ADR fixes the
  architecture), ADR-028 (adds `Bnormal_from_net_coil_currents` to `/problem`).
- **Context:** Phase 10 folds `regcoilPlot`, `compareRegcoil`, the coil cutters,
  and the Plotly port of `m20160811_01_plotCoilsFromRegcoil.m` into `regcoil`.
  Requirements: make **one plot type at a time** (not only the legacy fixed
  multi-panel dashboard) yet keep a **single "plot everything"** command; plot
  from a **live** run (no save) or a **saved** run with **no Fortran kernel and no
  expensive BLAS** (the ADR-028 promise); support the legacy plot menagerie —
  surface cross-sections, Pareto fronts (`f_B`/`f_K`, `max_Bnormal`/`max_K`, and
  the mixed pairs) for one run or several overlaid, the `(θ, ζ)` field maps
  (single-valued Φ, total Φ, K, Bnormal) across a λ scan, a 3D interactive view of
  plasma and/or winding surfaces, and a 3D interactive view of finite-thickness
  cut coils.
- **Decisions:**
  1. **Consume the object model only.** Every plotting function takes in-memory
     `PlasmaSurface` / `CoilSurface` / `Regcoil` / `Solution` / `SolutionScan`
     objects (or the `regcoil.load()` container), **never** file paths or raw
     `xarray`. Because `load()` reconstructs real objects via the ADR-028
     cheap-assembly split, a live object and a loaded object present one identical
     interface to the plot layer. File → object is the CLI's job, not the plot
     layer's.
  2. **Canonical free functions in `regcoil.plot`; thin delegating methods on the
     objects (hybrid).** `regcoil.plot.<fn>(obj, ...)` is the implementation;
     `obj.plot_*()` convenience methods delegate to it. matplotlib and plotly are
     imported **lazily inside** the plot functions/methods (the pattern already in
     `Surface.plot`), so `import regcoil` never imports a plotting library and the
     core solve keeps no plotting dependency.
  3. **Atomic functions return their axes/figure; dashboards are pure
     composition.** Every atomic function accepts `ax=` (matplotlib) or `fig=`
     (plotly) and returns it. The λ-scan panel grids and `plot.all()` create the
     layout and call the atomic functions into it — they contain **no** plotting
     logic of their own. This single invariant is what lets "one plot at a time,"
     the multi-panel grids, "overlay several runs on one axes," and "plot
     everything" all share code with zero duplication.
  4. **matplotlib for 2D, Plotly for 3D interactive.** 2D (cross-sections,
     Pareto/λ traces, `(θ, ζ)` contour maps) → matplotlib. 3D interactive
     (surfaces with a translucent-or-wireframe outer surface; finite-thickness cut
     coils) → Plotly (ADR-013). No backend-abstraction layer; each function
     documents which library's object it returns.
  5. **Default matplotlib figure size** is `regcoil.plot.DEFAULT_FIGSIZE =
     (14.5, 8.1)`, applied whenever a function must create its own figure (i.e.
     `ax=` was not supplied). Passing `ax=`/`fig=` overrides it entirely.
  6. **`cross_section` default angles.** `plot.cross_section(plasma, coil)`
     defaults to `phi = np.array([0, 0.5, 1, 1.5]) * np.pi / nfp` (half a field
     period, four slices — the legacy `regcoilPlot` set); a user may pass any
     `phi` array. Cross sections are computed by a new `Surface.cross_section(phi)
     -> (R, Z)` geometry method that derives `phi = atan2(y, x)` from the actual
     Cartesian `r` grid and interpolates each θ-line to the requested physical
     angle. This is correct for **both** `standard_toroidal_angle` values (ADR-025)
     — the plot layer never assumes constant-ζ == constant physical angle.
  7. **Composable 3D scene — `plot.plot_3d`.** A single 3D function plots **any
     subset** of {plasma surface, winding surface, cut coils} in one view via
     `plot.plot_3d(plasma=..., winding_surface=..., coils=...,
     winding_surface_style="wireframe")`; any argument may be omitted. The
     argument for the `CoilSurface` is `winding_surface=` (not `coil=`) so it does
     not collide with `coils=`, the discrete `CutCoils`. The render style applies
     specifically to the winding surface, so it is `winding_surface_style ∈
     {"wireframe", "translucent", "solid"}` (not a generic `outer=`).
  8. **Cutting is compute, separate from plotting.** `regcoil.cut(solution,
     coils_per_half_period, thickness=..., ...) -> CutCoils` contours the total Φ,
     maps contours to 3D curves via the coil Fourier modes, and optionally builds
     finite-thickness ribbons. `CutCoils.save_makegrid(path)` replaces
     `cutCoilsFromRegcoil`'s file write; `plot.plot_3d(coils=cut)` /
     `plot.coil_3d(cut)` render it. The irreversible file write stays out of the
     plot path.
  9. **A `Bnormal_from_net_coil_currents` panel needs no kernel on load** because
     ADR-028 (amended) stores that λ-independent plasma-grid array in `/problem`.
     Without this, that one legacy panel would have forced an operator rebuild on a
     loaded run, violating the ADR-028 promise.
  10. **One `regcoil` CLI with subcommands** (`regcoil plot|compare|cut|...`),
      each a thin `load() → plot function(s) → show/save` wrapper. No separate
      `regcoilPlot`/`compareRegcoil` console names. See [API.md](API.md#command-line-interface).
- **Consequences:**
  - New `src/regcoil/plot.py` (or `plot/` package) holding the free functions and
    `DEFAULT_FIGSIZE`; new `src/regcoil/cut.py` with `cut()` / `CutCoils`; a
    `Surface.cross_section` method; convenience `plot_*` methods on the objects.
  - `Surface.plot` (the current minimal matplotlib 3D stub) is superseded by the
    Plotly `plot_3d` path; the matplotlib stub may stay as a lightweight fallback
    or be removed.
  - Plotting smoke tests (Phase 11) exercise the functions on non-interactive
    backends against a small live run and a round-tripped saved run.
  - `plotly` is already an allowed runtime dependency (ADR-011); no new dep.

---

## ADR-030: Documentation stack details (Phase 12)

- **Date:** 2026-07-20
- **Status:** accepted
- **Supersedes:** ADR-014 (keeps its three firm decisions — publish on Read the
  Docs, retire `manual/` as canonical, delete `publish_manual.yml` outright — and
  resolves the choices it left open).
- **Context:** Phase 12 replaces the LaTeX `manual/` with a Read the Docs manual.
  Requirements gathered for this phase: a short quickstart in `README.md`; pages
  for installation, typical usage, **every** plot type, and an API reference; an
  API-reference presentation in the style of yancc's `api.html` (grouped
  `autosummary` tables with per-entry one-line descriptions and links to source);
  and a hard guarantee that **code shown in the docs actually runs as written**,
  ideally by executing the doc pages the way DESC executes notebooks in its build.
  DESC's own model (checked-in `.ipynb` with committed outputs and
  `nbsphinx_allow_errors = True`) meets the "notebooks in the docs" goal but not
  the "proven to run" goal, and its multi-MB JSON blobs diff and merge badly. A
  pyREGCOIL-specific wrinkle: constructing and *solving* a `Regcoil` calls the
  compiled Fortran kernel, but by ADR-028's cheap-assembly split, `load()`-ing a
  saved run and plotting it needs **no** Fortran kernel and no expensive BLAS — so
  most doc pages can execute with zero build toolchain if they load a shipped run.
- **Decisions:**
  1. **Sphinx, not MkDocs.** The requested API look is a Sphinx `autosummary`
     page, and the surrounding stellarator-code ecosystem (yancc, DESC, SIMSOPT)
     is uniformly Sphinx, so contributors and users are already fluent in it.
  2. **Theme: `furo`.** Modern look, first-class dark mode, pairs cleanly with
     MyST-NB. (Chosen over `sphinx-book-theme`.)
  3. **API reference via `sphinx.ext.autosummary` with `:recursive:` +
     `:toctree: _api/`, grouped by topic, one line per entry — the yancc pattern.**
     Groups: Surfaces (`Surface`, `FourierSurface`, `PlasmaSurface`,
     `CoilSurface`), Problem & solution (`Regcoil`, `Solution`, `SolutionScan`),
     Plotting (`regcoil.plot.*`), Coil cutting (`regcoil.cut`, `CutCoils`),
     Save/load (`regcoil.save`/`load`). Docstrings are NumPy-style, rendered by
     `sphinx.ext.napoleon`; `sphinx.ext.linkcode` (or `sphinx_github_style`, as in
     yancc) adds a per-object **[source]** link to the exact lines on GitHub.
  4. **Executed pages via MyST-NB, failing the build on any cell error.**
     `nb_execution_mode = "cache"` (jupyter-cache, so unchanged pages are not
     re-run) and `nb_execution_raise_on_error = True` (a broken cell fails the
     build — this is what makes "the code runs" a checked invariant, unlike DESC's
     `allow_errors`). This is the primary runnable-code mechanism, replacing
     `nbsphinx`.
  5. **Doc-page source is jupytext-paired MyST-Markdown (`.md`), with no committed
     cell outputs.** Plain-text, reviewable diffs; still openable/editable as a
     notebook in Jupyter via the jupytext pairing; outputs are produced at build
     time by decision 4, never stored in git. This deliberately avoids DESC's
     multi-MB `.ipynb`-with-outputs blobs.
  6. **A second safety net in CI: `pytest --nbmake` over the doc pages, on the
     normal test matrix**, plus a `sphinx-build -W` (warnings-as-errors) smoke
     build. This catches doc breakage in ordinary PRs across Python versions
     before Read the Docs ever builds, and does not depend on the RTD service.
  7. **The `README.md` quickstart snippet is validated by `sphinx.ext.doctest`**
     (`make -C docs doctest` in CI). The README carries a short quickstart; a
     fuller, executed quickstart lives in `docs/quickstart.md`. `doctest` is used
     only for that short snippet — tutorials use MyST-NB (decision 4), not doctest.
  8. **Read the Docs compiles the Fortran extension** so executed pages run the
     real code path: `.readthedocs.yaml` sets `build.apt_packages: [gfortran,
     libopenblas-dev]` and pip-installs the package with its docs extra. RTD is the
     canonical published build; `sphinx.fail_on_warning: true`. (The alternative —
     execute in GitHub Actions and have RTD only render pre-executed artifacts — is
     the documented fallback if RTD build time becomes a problem; it is not the
     default because it reintroduces stored-artifact plumbing.)
  9. **Ship one or two small saved runs (`.nc`) under `docs/`; kernel-free pages
     `load()` + plot them.** By the ADR-028 promise, those pages (installation,
     the plotting gallery, most of "typical usage") execute with no Fortran and no
     BLAS. Only a single "construct + solve from scratch" tutorial requires the
     compiled extension, so decision 8's toolchain is exercised on a minimal
     surface and the bulk of the docs stay cheap to build.
- **Consequences:**
  - New `docs/` layout (see [PHASES.md](PHASES.md#phase-12) Phase 12): `conf.py`,
    `index.md`, `installation.md`, `quickstart.md`, `usage.md`, `plotting.md`,
    `api.md`, `theory.md`, `_static/`, and small `*.nc` sample run(s).
  - New `.readthedocs.yaml` (apt gfortran + libopenblas, docs extra, fail on
    warning); a `docs` extra in `pyproject.toml`'s optional dependencies:
    `sphinx`, `furo`, `myst-nb`, `jupytext`, `sphinx-copybutton`,
    `sphinx-github-style` (or `linkcode`), `nbmake` (dev/CI). Docs deps are
    docs/dev-only, never runtime (carries ADR-014's separation).
  - CI gains a docs job: `pytest --nbmake docs`, `sphinx-build -W`, and
    `make -C docs doctest`.
  - `manual/` is retired as canonical and `.github/workflows/publish_manual.yml`
    is deleted (ADR-014 / Phase 12 unchanged); the theory page ports the essential
    physics from `manual/*.tex`.
  - Untracked root scratch scripts (`run_example.py`, `run_fast_example.py`,
    `run_slow_example.py`, `test_cross_section_plot.py`) are candidate raw material
    for the usage/plotting pages and should be removed from the repo root once
    their content is folded into docs (or into `examples/`).
