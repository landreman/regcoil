### REGCOIL

A regularized current potential method for fast computation of the shapes of stellarator coils

![](https://github.com/landreman/regcoil/blob/master/docs/m20170111_01_compareNescoilToRegcoilCoils.png)

This program is described in the paper
`M Landreman, "An improved current potential method for fast computation of stellarator coil shapes," Nuclear Fusion 57, 046003 (2017)`,
which is available in this repository and also [at arXiv:1609.04378](https://arxiv.org/pdf/1609.04378.pdf).

This repository provides a pip-installable python package, with compiled kernels
for the computationally intensive steps.
The original fortran implementation of the REGCOIL algorithm is available in the
companion repository
https://github.com/landreman/regcoil_fortran

### Installation

This package is available on pypi (which provides pre-compiled wheels), so it
can be installed in the usual way with pip:

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
pip install "meson-python>=0.16.0"
pip install .[dev]
```

For editable installs from source, it is necessary to include the ``--no-build-isolation``
flag due to a quirk of the meson build system:

```bash
pip install --no-build-isolation -e .[dev]
```

From the cloned repository, tests can be run using

```bash
pytest
```

### Quickstart

<!-- quickstart-start -->

```python
import regcoil

ds = regcoil.examples("W7-X")  # Or "NCSX"
# ds then provides paths to a vmec wout file, simsopt virtual casing file,
# stellopt bnorm file, and coil winding surface in nescin format.

# Define the plasma boundary surface:
plasma = regcoil.PlasmaSurface.from_wout(ds.wout, ntheta=64, nzeta=64)
# Assign B_normal data associated with the plasma current:
plasma.set_bnormal_from_virtual_casing(ds.vcasing)

# Define a coil winding surface:
coil = regcoil.CoilSurface.from_uniform_offset(
    plasma, separation=0.3, ntheta=64, nzeta=64, mpol=12, ntor=12
)

problem = regcoil.Regcoil(plasma, coil, mpol_potential=12, ntor_potential=12)
solution = problem.solve(lam=1e-14)
# You can also solve for a range of lambda values, or search for lambda that 
# results in a desired f_B, f_K, max_K, etc
print(f"f_B = {solution.f_B:.1e}, f_K = {solution.f_K:.1e}")
# f_B = 2.7e-01, f_K = 1.1e+15

# Results can then be saved and plotted.
```

<!-- quickstart-end -->

Full documentation, covering installation, typical usage, saving/loading data, plotting, and
the API reference, are available at [regcoil.readthedocs.io](https://regcoil.readthedocs.io).

