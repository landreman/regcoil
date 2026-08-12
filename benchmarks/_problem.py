"""The problem the whole benchmark suite is built on, plus its size.

Every benchmark works on the bundled W7-X example (the README quickstart case),
so the suite exercises one physically meaningful configuration end to end: read
the VMEC boundary, build a winding surface, assemble the operators, solve, cut
coils, save/load.

The resolution is deliberately modest.  The benchmarks run under CodSpeed's
CPU-simulation instrument, which executes the measured code inside Valgrind at
roughly two orders of magnitude of slowdown, so the useful size is the smallest
one that still spends its time where the production-size problem does: the
Fortran kernel, the BLAS level-3 calls, the LAPACK eigensolve, and the numpy
geometry.

Changing any of the constants here invalidates the historical CodSpeed series
for these benchmarks, so change them deliberately.
"""

from __future__ import annotations

import numpy as np

import regcoil

#: Grid used for both the plasma and the coil surface.
NTHETA = 32
NZETA = 32
#: Fourier resolution of the winding surface.
MPOL = 8
NTOR = 8
#: Fourier resolution of the current potential (144 basis functions).
MPOL_POTENTIAL = 8
NTOR_POTENTIAL = 8
#: Plasma-to-coil-surface separation, in meters.
SEPARATION = 0.3
#: Regularization weight for the single-lambda solves.
LAM = 1e-14
#: Lambda values for the scan benchmark.
LAMBDAS = np.logspace(-19, -11, 32)
#: Coils cut per half field period.
COILS_PER_HALF_PERIOD = 5
#: Ribbon width, in meters, for the finite-thickness cut.
THICKNESS = 0.05

#: Coarser grid for `standard_toroidal_angle=True`, whose per-grid-point root
#: solve costs several times more than the along-normal construction.
NTHETA_COARSE = 16
NZETA_COARSE = 16


def make_plasma(w7x):
    """The W7-X plasma boundary with virtual-casing B_normal."""
    surface = regcoil.PlasmaSurface.from_wout(w7x.wout, ntheta=NTHETA, nzeta=NZETA)
    surface.set_bnormal_from_virtual_casing(w7x.vcasing)
    return surface


def make_coil(plasma):
    """The uniform-offset winding surface used by the solver benchmarks."""
    return regcoil.CoilSurface.from_uniform_offset(
        plasma, separation=SEPARATION, ntheta=NTHETA, nzeta=NZETA,
        mpol=MPOL, ntor=NTOR,
    )


def make_problem(plasma, coil):
    """Assemble the operators: Fortran kernel, `matrix_B`/`matrix_K`, eigensolve."""
    return regcoil.Regcoil(
        plasma, coil,
        mpol_potential=MPOL_POTENTIAL, ntor_potential=NTOR_POTENTIAL,
    )
