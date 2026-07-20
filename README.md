### REGCOIL

A regularized current potential method for fast computation of the shapes of stellarator coils

![](https://github.com/landreman/regcoil/blob/master/manual/m20170111_01_compareNescoilToRegcoilCoils.png)

This program is described in the paper
`M Landreman, "An improved current potential method for fast computation of stellarator coil shapes," Nuclear Fusion 57, 046003 (2017)`,
which is available in this repository and also [at arXiv:1609.04378](https://arxiv.org/pdf/1609.04378.pdf).
For further documentation, see the user manual [here](http://landreman.github.io/regcoil/regcoilManual.pdf).

Python packaging / migration plans: [`docs/migration/`](docs/migration/OVERVIEW.md).

### Install (Python package)

Requires **gfortran**, **BLAS**, and **OpenMP** (the Python package's stateless Fortran kernels need no LAPACK and no NetCDF, as of Phase 9). The legacy executable below still needs **LAPACK** and **NetCDF C + Fortran** until Phase 13. On macOS Homebrew Apple Silicon see [`docs/migration/LOCAL.md`](docs/migration/LOCAL.md).

Editable installs with meson-python need build tools in the *same* environment (pip’s isolated build env is discarded after install):

```bash
pip install -e ".[dev]" --no-build-isolation
```

If build tools are not yet installed:

```bash
pip install meson ninja "meson-python>=0.16" pytest
pip install -e ".[dev]" --no-build-isolation
```

Non-editable install (fine for CI / one-shot use; uses build isolation):

```bash
pip install ".[dev]"
```

Run tests (requires the installed python package):

```bash
pytest
```

### Build and test (legacy executable)

The Fortran sources live under `fortran/`. The root `makefile` still builds the `regcoil` binary at the repo root for example regressions:

```bash
make
make test
```

On macOS with Homebrew NetCDF (Apple Silicon), the default makefile paths under `/usr/local` are often wrong—see [`docs/migration/LOCAL.md`](docs/migration/LOCAL.md) for the `EXTRA_*` flag overrides that work locally.
