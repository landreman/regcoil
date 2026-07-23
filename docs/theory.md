# Theory

This page summarizes the physical model and notation used throughout the
API. For more details, see
[Landreman, *Nuclear Fusion* 57, 046003 (2017)](https://arxiv.org/pdf/1609.04378.pdf)
(also included in the repository as a PDF).

## The current potential

REGCOIL represents the coil currents as a sheet current on a toroidal
"winding" surface, described by a scalar current potential $\Phi(\theta,
\zeta)$ (poloidal angle $\theta$, toroidal angle $\zeta$). The potential
splits into a secular part carrying the net poloidal and toroidal currents,
plus a single-valued, doubly-periodic part represented by a finite Fourier
series:

$$
\Phi(\theta, \zeta) = \frac{G}{2\pi}\zeta + \frac{I}{2\pi}\theta +
\Phi_{sv}(\theta, \zeta),
$$

where $G$ is the net poloidal current (`net_poloidal_current`) and $I$ is the
net toroidal current (`net_toroidal_current`), both in Amperes. The
single-valued part $\Phi_{sv}$ is expanded in the Fourier modes selected by
`mpol_potential`/`ntor_potential`; its coefficients
(`single_valued_current_potential_mn`) are exactly what `Regcoil.solve`
solves for.

The surface current density follows from the potential by

$$
\mathbf{K} = \hat{\mathbf{n}} \times \nabla \Phi,
$$

available as `Solution.current_density()`.

## The regularized least-squares problem

Given the plasma surface's target normal field $B_{n}^{\mathrm{target}}$
(typically from in-plasma current, loaded via
`PlasmaSurface.set_bnormal_from_bnorm_file` or
`PlasmaSurface.set_bnormal_from_virtual_casing`), the coils must produce a normal
field that cancels it. REGCOIL minimizes a weighted sum of two objectives
over the Fourier coefficients of $\Phi_{sv}$, for a chosen regularization
parameter $\lambda$:

$$
f = f_B + \lambda\, f_K,
\qquad
f_B = \int d^2a\, \big(\mathbf{B}\cdot\hat{\mathbf{n}}\big)^2,
\qquad
f_K = \int d^2a\, K^2,
$$

where the $f_B$ integral is over the plasma surface and the $f_K$
integral is over the coil surface. Because both terms are quadratic in the
unknown Fourier coefficients, the minimizer is a linear solve; scanning many
values of $\lambda$ (`Regcoil.scan`) or bisecting for a target
(`Regcoil.solve_for_target`) reuses one eigendecomposition of the problem
(`scipy.linalg.eigh`) rather than repeating an expensive assembly step (see
the `Regcoil` API reference).

In the Python API these two objectives are named `Solution.f_B` and
`Solution.f_K` (renamed from the paper's/legacy code's `chi2_B`/`chi2_K`). `lam=0`
recovers the least-squares fit with no regularization; `lam=np.inf` recovers
the (well-defined) heavily-regularized limit.

## The trade-off

Smaller $\lambda$ drives $f_B \to 0$ (better field accuracy) at the cost
of larger, more convoluted currents (larger $f_K$, `max_K`, `rms_K`);
larger $\lambda$ favors simpler, lower-current coils at the cost of field
accuracy. Plotting `f_B` against `f_K` across a scan traces out this
Pareto front — see [`regcoil.plot.pareto`](plotting.md#pareto-fronts).

## From continuous current to discrete coils

A solved current potential is continuous; `regcoil.cut.cut` turns it into a
finite number of discrete coils by contouring $\Phi$ at evenly-spaced levels,
each contour becoming one coil that carries the current between adjacent
levels. See [Cutting into discrete coils](plotting.md#cutting-into-discrete-coils).

## Units

All inputs and outputs are SI: meters, Tesla, Amperes, and combinations thereof.

## Parallelization

The Fortran kernels (matrix assembly) are parallelized with OpenMP; matrix
products and the eigendecomposition use whatever BLAS/LAPACK NumPy/SciPy are
linked against. There is no MPI; REGCOIL runs on a single node, using however
many threads `OMP_NUM_THREADS` (and your BLAS library's own thread control)
allow.
