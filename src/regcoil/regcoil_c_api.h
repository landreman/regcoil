#ifndef REGCOIL_C_API_H
#define REGCOIL_C_API_H

#ifdef __cplusplus
extern "C" {
#endif

void regcoil_c_set_verbose(int flag);
int regcoil_c_nlambda(void);
int regcoil_c_setup(const char *path);
int regcoil_c_build_matrices(void);
int regcoil_c_prepare_solve(void);
int regcoil_c_solve_ilambda(int ilambda, double *chi2_b, double *chi2_k,
                            double *max_bn, double *max_k);
int regcoil_c_solve_lambda(double lam, double *chi2_b, double *chi2_k,
                           double *max_bn, double *max_k);

#ifdef __cplusplus
}
#endif

#endif /* REGCOIL_C_API_H */
