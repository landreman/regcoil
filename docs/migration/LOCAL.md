# Local build and test notes

Practical notes for developers and agents running the **current** Fortran executable path (`make` / `make test`) and the **pip-installed** Python package (`regcoil._core`). Packaging/CI supersede the makefile for the supported import path (Phase 4+).

Fortran sources: `fortran/` (including `fortran/mini_libstell/`). The root `makefile` uses `VPATH` and still places the `regcoil` binary at the repo root.

## Dependencies

- **Fortran:** `gfortran` (or Open MPI’s `mpif90` wrapper), BLAS/LAPACK, OpenMP
- **NetCDF:** C + Fortran libraries (`netcdf`, `netcdf-fortran`) until Phase 8
- **Python (examples):** `numpy`, `scipy` (tests use `scipy.io.netcdf_file`)
- **Python (package):** see root `pyproject.toml` / `pip install ".[dev]"` (meson-python; builds `regcoil._core`)

On this machine’s default Cursor conda env, use **`20231204-02-desc`** for Python (see workspace rules).

## local environment

For running python commands and tests on my macbook, you may use the conda env 20260711-01-opt. For installing this regcoil package or other packages, use a virtual environment to avoid messing up the env 20260711-01-opt.

**macOS + Homebrew NetCDF:** create the venv with Homebrew (or another non-conda) Python, e.g. `/opt/homebrew/bin/python3 -m venv .venv`. A venv based on a conda env often loads conda’s HDF5 first and then `import regcoil._core` fails with a missing `H5T_IEEE_F16BE_g` symbol when the extension was linked against Homebrew NetCDF/HDF5.

## macOS (Homebrew / Apple Silicon)

The makefile’s default host stanza still points at `/usr/local`. On Homebrew Apple Silicon that is wrong; NetCDF lives under `/opt/homebrew`.

**Important:** the makefile *assigns* `EXTRA_COMPILE_FLAGS` / `EXTRA_LINK_FLAGS` in host `ifeq` blocks, so shell `export` of those variables is ignored. Pass them on the **`make` command line**:

```bash
make -j4 \
  EXTRA_COMPILE_FLAGS="-fopenmp -I/opt/homebrew/include -ffree-line-length-none -O0 -g -fallow-argument-mismatch" \
  EXTRA_LINK_FLAGS="-fopenmp -L/opt/homebrew/lib -lnetcdff -lnetcdf -framework Accelerate"

make test \
  EXTRA_COMPILE_FLAGS="-fopenmp -I/opt/homebrew/include -ffree-line-length-none -O0 -g -fallow-argument-mismatch" \
  EXTRA_LINK_FLAGS="-fopenmp -L/opt/homebrew/lib -lnetcdff -lnetcdf -framework Accelerate"
```

Confirm NetCDF Fortran is installed (`netcdf.inc` should exist), e.g.:

```bash
ls /opt/homebrew/include/netcdf.inc
# or: brew list netcdf-fortran
```

For **pip** builds on Apple Silicon, ensure Homebrew `pkg-config` can see NetCDF (`brew install pkg-config netcdf netcdf-fortran`). Meson uses `dependency('netcdf')` / `netcdf-fortran` and Accelerate when LAPACK pkg-config is missing.

## Tests today

- `make test` → `examples/runExamples.py` (discovers `examples/*/tests.py` + matching `regcoil_in.*`)
- Expect several minutes for the full suite; success ends with `ALL TESTS THAT WERE RUN WERE PASSED SUCCESSFULLY` and `numFailures: 0`
- Package tests: create/use a **venv**, then
  `pip install ".[dev]"` (or editable with `--no-build-isolation`) and `pytest`.
  Includes `import regcoil` / `RegcoilProblem`, a one-λ axisymmetry parity test, and a two-instance non-interference test (`tests/unit/test_core_one_lambda.py`).

## GitHub Actions (Phase 3–4 / ADR-016, ADR-017)

Workflow: `.github/workflows/ci.yml`. Builds the legacy executable **and** installs the Python package via **pip** (compiling `regcoil._core`), then runs pytest on `ubuntu-latest` and `macos-latest`. Full example regressions are **not** in CI yet.

Reproduce the CI makefile hosts locally:

```bash
# Linux (needs gfortran + BLAS/LAPACK + NetCDF Fortran)
make -j4 MACHINE=github_ubuntu

# macOS Homebrew
brew install gcc netcdf netcdf-fortran
make -j4 MACHINE=github_macos
```

## Other hosts

Use the makefile’s named host stanzas (`MACHINE=…`, `CLUSTER=…`, NERSC, etc.) when they match your machine; see comments at the top of `makefile`.
