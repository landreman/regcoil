"""Golden reference values for the stateless Fortran kernels
(`regcoil._core.build_g_and_h`, `.build_inductance`, `.uniform_offset_surface`),
generated from the legacy compiled Fortran and hardcoded here.

g/h case: a manufactured (non-physical) small grid -- ntheta_plasma=3,
nzeta_plasma=4, ntheta_coil=5, nzeta_coil=6, nfp=2 -- with deterministic
sin/cos geometry (see `test_kernels.py::_manufactured_geometry`) fed directly
into the legacy `regcoil_build_matrices` (mpol_potential=1, ntor_potential=0,
symmetry_option=1, so the single basis function is sin(theta_coil)) and
compared digit-for-digit against `regcoil._core.build_g_and_h`.

Offset-surface case: a plasma shape with a helical (m=1, n=1) bump --
xm=[0,1,1], xn=[0,0,3] (xn already includes the nfp=3 factor), rmnc=[5, 1,
0.15], zmns=[0, 1, 0.1] -- run through the legacy `regcoil_init_coil_surface`
(geometry_option_coil=2, separation=0.4, ntheta_coil=6, nzeta_coil=5,
max_mpol_coil=3, max_ntor_coil=2, no filtering) and compared against
`regcoil._core.uniform_offset_surface`. The bump forces a genuine root solve
in `regcoil_fzero` (a circular cross section needs none), so this also
exercises regcoil_fzero.f.
"""

import numpy as np

# --- build_g_and_h / build_inductance golden case ---------------------------

G_GOLDEN = np.array([
    9.0277637934037047e-10, 9.0997313026554296e-10, 8.2709276396528163e-10,
    7.6262079247585073e-10, 8.3479092139693841e-10, 8.3556855413484085e-10,
    6.1879482831302066e-10, 7.4762479378537296e-10, 8.3580967433387366e-10,
    5.1344050885466239e-10, 6.6873181361786127e-10, 8.1949873480150019e-10,
])

H_GOLDEN = np.array([
    1.1132227973270049e-11, 8.4073577988661823e-12, 3.2454912900852670e-12,
    7.5793754387724669e-12, 5.8563808601848271e-12, 2.0815329247749119e-12,
    4.6053328924349743e-12, 4.5267096166502054e-12, 2.7642137673236110e-12,
    2.4837051847633119e-12, 4.1408014808101436e-12, 4.8732931516000988e-12,
])

# --- uniform_offset_surface golden case --------------------------------------

XM_OUT = np.array([0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3])
XN_OUT = np.array([0, 3, 6, -6, -3, 0, 3, 6, -6, -3, 0, 3, 6, -6, -3, 0, 3, 6])

RMNC_OUT = np.array([
    4.9997798173568251e+00, 1.0860320125988424e-05, 1.3142314680417781e-04,
    2.9920567272547043e-03, -2.4895187634708830e-02, 1.3997223183639509e+00,
    1.7006474521150661e-01, -1.9998626159365997e-03, -1.8997200949985477e-05,
    -2.3295683872872478e-05, -2.1952286863553832e-04, -4.1313268765804225e-05,
    3.8102819847041838e-04, -2.2804497793529312e-05, 2.5787721653984692e-03,
    -9.9600538727600088e-04, 2.5787721653985525e-03, -2.2804497793364628e-05,
])
RMNS_OUT = np.zeros(18)
ZMNC_OUT = np.zeros(18)
ZMNS_OUT = np.array([
    0.0000000000000000e+00, 7.5554790417764581e-05, 4.4144138780256647e-04,
    3.0078292748591537e-03, -2.4680004690783894e-02, 1.3983131455402795e+00,
    1.3006920830554280e-01, -2.6419942668672795e-03, -1.7434314236400742e-04,
    -4.7186117656571746e-05, -8.2986077214094095e-05, 1.0116768359479231e-05,
    3.5964241191324269e-04, -6.0345934099814614e-05, -2.3683869490397085e-03,
    -4.9425075307514537e-17, 2.3683869490397423e-03, 6.0345934099756331e-05,
])
