# Installation

`regcoil` is a Python package with a small compiled Fortran extension
(`regcoil._core`) built by [meson-python](https://mesonbuild.com/meson-python).

## Install with pip

This package is available on pypi (which provides pre-compiled wheels), so it can be installed in the usual way with pip:
```bash
pip install regcoil
```

If you want updates that are more recent than the most recent release, or if you want to edit the
source code,
you can also install regcoil from source.  To do so, first clone the repository and then
pip-install from the local repository:

```bash
git clone https://github.com/landreman/regcoil.git
cd regcoil
pip install .[dev]
```

For editable installs from source, it is necessary to include the ``--no-build-isolation``
flag due to a quirk of the meson build system:

```bash
pip install --no-build-isolation -e .[dev]
```

To verify the extension built and imports correctly:

```bash
python -c "import regcoil._core; print('regcoil._core OK')"
```

## Running the test suite

The tests can be run using

```bash
pytest
```

Several high-resolution (`ntheta_plasma=nzeta_plasma=128`) regression cases are marked `slow`
and can be skipped if desired:

```bash
pytest -m "not slow"
```
