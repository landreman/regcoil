# Installation

`regcoil` is a Python package with a small compiled Fortran extension
(`regcoil._core`) built by [meson-python](https://mesonbuild.com/meson-python).

## Requirements

Building the extension needs:

- gfortran
- BLAS (on macOS, the Accelerate framework is used automatically if no
  BLAS `pkg-config` file is found)
- OpenMP

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

The high-resolution (`ntheta_plasma=nzeta_plasma=128`) regression cases are marked `slow`
and can be skipped if desired:

```bash
pytest -m "not slow"
```
