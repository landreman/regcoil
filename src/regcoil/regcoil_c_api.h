#ifndef REGCOIL_C_API_H
#define REGCOIL_C_API_H

#ifdef __cplusplus
extern "C" {
#endif

/* Compiled stateless kernels. All arrays are caller-allocated,
 * Fortran-contiguous (column-major), float64 (int32 for mode-number
 * arrays), and passed as raw pointers; a nonzero return value is an error
 * code (see the corresponding Fortran `info`). No opaque handle, no
 * persistent state: concurrent calls with different sizes do not
 * interfere. */

int regcoil_c_omp_max_threads(void);

int regcoil_c_build_inductance(
    int ntheta_plasma, int nzeta_plasma, int ntheta_coil, int nzeta_coil, int nfp,
    const double *r_plasma, const double *normal_plasma,
    const double *r_coil, const double *normal_coil,
    const double *drdtheta_coil, const double *drdzeta_coil,
    double net_poloidal_current, double net_toroidal_current,
    double dtheta_coil, double dzeta_coil,
    double *inductance, double *h);

int regcoil_c_build_g_and_h(
    int ntheta_plasma, int nzeta_plasma, int ntheta_coil, int nzeta_coil, int nfp, int nbf,
    const double *r_plasma, const double *normal_plasma,
    const double *r_coil, const double *normal_coil,
    const double *drdtheta_coil, const double *drdzeta_coil,
    const double *basis_functions,
    double net_poloidal_current, double net_toroidal_current,
    double dtheta_coil, double dzeta_coil,
    double *g, double *h);


#ifdef __cplusplus
}
#endif

#endif /* REGCOIL_C_API_H */
