# Installation

`regcoil` is a Python package with a small compiled Fortran extension
(`regcoil._core`) built by [meson-python](https://mesonbuild.com/meson-python).

## Requirements

Building the extension needs:

- **gfortran**
- **BLAS** (on macOS, the Accelerate framework is used automatically if no
  BLAS `pkg-config` file is found)
- **OpenMP**

No LAPACK and no NetCDF are required to build or import the package: the
stateless Fortran kernels only need BLAS, and all file I/O (VMEC `wout`,
BNORM, NESCIN, and `regcoil`'s own save/load format) is pure Python.

## Install with pip

Editable installs with meson-python need the build tools (`meson`, `ninja`,
`meson-python`) installed in the *same* environment, since pip's isolated
build environment is discarded after the install finishes:

```bash
pip install meson ninja "meson-python>=0.16" pytest
pip install -e ".[dev]" --no-build-isolation
```

A plain, non-editable install (fine for CI or one-off use) can rely on pip's
build isolation instead, and does not need the tools above pre-installed:

```bash
pip install "regcoil[dev]"
```

Verify the extension built and imports correctly:

```bash
python -c "import regcoil._core; print('regcoil._core OK')"
```

## Running the test suite

```bash
pytest
```

The high-resolution (`ntheta_plasma=128`) regression cases are marked `slow`
and skipped by default:

```bash
pytest -m "not slow"   # matches CI on Linux
pytest                 # run everything, including slow cases
```

## macOS (Homebrew, Apple Silicon)

See [`docs/migration/LOCAL.md`](https://github.com/landreman/pyREGCOIL/blob/master/docs/migration/LOCAL.md)
in the repository for Homebrew-specific notes, including the venv/conda
interaction that can otherwise cause `import regcoil._core` to fail with a
missing HDF5 symbol.

## The legacy Fortran executable

The repository also still builds a standalone Fortran executable via the root
`makefile`, used for regression comparisons during the ongoing Python
migration ([`docs/migration/`](https://github.com/landreman/pyREGCOIL/tree/master/docs/migration)).
It requires LAPACK and NetCDF (C + Fortran) in addition to the requirements
above:

```bash
make
make test
```

This executable is not part of the pip-installed `regcoil` package and is not
required to use `regcoil` from Python.
