# Local build and test notes

Practical notes for developers and agents running the **current** Fortran executable path (`make` / `make test`). Packaging/CI will supersede parts of this in later phases.

## Dependencies

- **Fortran:** `gfortran` (or Open MPI’s `mpif90` wrapper), BLAS/LAPACK, OpenMP
- **NetCDF:** C + Fortran libraries (`netcdf`, `netcdf-fortran`) until Phase 8
- **Python (examples):** `numpy`, `scipy` (tests use `scipy.io.netcdf_file`)

On this machine’s default Cursor conda env, use **`20231204-02-desc`** for Python (see workspace rules).

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

## Tests today

- `make test` → `examples/runExamples.py` (discovers `examples/*/tests.py` + matching `regcoil_in.*`)
- Expect several minutes for the full suite; success ends with `ALL TESTS THAT WERE RUN WERE PASSED SUCCESSFULLY` and `numFailures: 0`

## Other hosts

Use the makefile’s named host stanzas (`MACHINE=…`, `CLUSTER=…`, NERSC, etc.) when they match your machine; see comments at the top of `makefile`.
