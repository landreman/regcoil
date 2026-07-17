### REGCOIL

A regularized current potential method for fast computation of the shapes of stellarator coils

![](https://github.com/landreman/regcoil/blob/master/manual/m20170111_01_compareNescoilToRegcoilCoils.png)

This program is described in the paper
`M Landreman, "An improved current potential method for fast computation of stellarator coil shapes," Nuclear Fusion 57, 046003 (2017)`,
which is available in this repository and also [at arXiv:1609.04378](https://arxiv.org/pdf/1609.04378.pdf).
For further documentation, see the user manual [here](http://landreman.github.io/regcoil/regcoilManual.pdf).

Python packaging / migration plans: [`docs/migration/`](docs/migration/OVERVIEW.md).

### Install (Python package scaffold)

Requires a Fortran compiler only once the extension is wired (Phase 4). For now the install is pure Python.

Editable installs with meson-python need build tools in the *same* environment (pip’s isolated build env is discarded after install):

```bash
pip install -e ".[dev]" --no-build-isolation
```

If build tools are not yet installed:

```bash
pip install meson ninja "meson-python>=0.16" pytest
pip install -e ".[dev]" --no-build-isolation
```

Non-editable install (fine for CI / one-shot use):

```bash
pip install .
# or with optional NetCDF4 (ADR-004 still open):
pip install ".[netcdf]"
```

Smoke test:

```bash
pytest
```

Allowed runtime dependencies: `numpy`, `scipy`, `matplotlib`, `f90nml`, `plotly` (optional `netCDF4`). No SIMSOPT/DESC.

Build backend: **meson-python** (ADR-002). Layout: `src/regcoil/`, `fortran/`, `tests/`, `examples/`, `docs/`.

### Build and test (legacy executable)

The Fortran sources live under `fortran/`. The root `makefile` still builds the `regcoil` binary at the repo root for example regressions:

```bash
make
make test
```

On macOS with Homebrew NetCDF (Apple Silicon), the default makefile paths under `/usr/local` are often wrong—see [`docs/migration/LOCAL.md`](docs/migration/LOCAL.md) for the `EXTRA_*` flag overrides that work locally.
