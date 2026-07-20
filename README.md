### REGCOIL

A regularized current potential method for fast computation of the shapes of stellarator coils

![](https://github.com/landreman/regcoil/blob/master/manual/m20170111_01_compareNescoilToRegcoilCoils.png)

This program is described in the paper
`M Landreman, "An improved current potential method for fast computation of stellarator coil shapes," Nuclear Fusion 57, 046003 (2017)`,
which is available in this repository and also [at arXiv:1609.04378](https://arxiv.org/pdf/1609.04378.pdf).

**Full documentation, including installation, a tutorial, every plot type, and
the API reference, is at [regcoil.readthedocs.io](https://regcoil.readthedocs.io).**

Python packaging / migration plans: [`docs/migration/`](docs/migration/OVERVIEW.md).

### Quickstart

<!-- quickstart-start -->

```python
import regcoil

ds = regcoil.examples("NCSX")
plasma = regcoil.PlasmaSurface.from_wout(ds.wout, ntheta=24, nzeta=24)
plasma.set_bnormal_from_bnorm_file(ds.bnorm)
coil = regcoil.CoilSurface.from_uniform_offset(
    plasma, separation=0.3, ntheta=24, nzeta=24, mpol=8, ntor=8
)

problem = regcoil.Regcoil(plasma, coil, mpol_potential=8, ntor_potential=8)
solution = problem.solve(lam=1e-14)
print(f"f_B = {solution.f_B:.1e}, f_K = {solution.f_K:.1e}")
# f_B = 1.9e-01, f_K = 9.7e+13
```

<!-- quickstart-end -->

`regcoil.examples` bundles a handful of ready-to-use VMEC/BNORM/NESCIN datasets
(`"NCSX"`, `"W7-X"`, ...) so the snippet above runs with no external files. See
[the quickstart page](https://regcoil.readthedocs.io/en/latest/quickstart.html)
for a fuller walkthrough (λ scans, plotting, saving/loading, cutting coils).

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
