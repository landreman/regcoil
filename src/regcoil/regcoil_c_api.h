#ifndef REGCOIL_C_API_H
#define REGCOIL_C_API_H

#ifdef __cplusplus
extern "C" {
#endif

typedef struct regcoil_problem regcoil_problem; /* opaque */

regcoil_problem *regcoil_c_create(void);
void regcoil_c_destroy(regcoil_problem *handle);
void regcoil_c_set_verbose(regcoil_problem *handle, int flag);
int regcoil_c_nlambda(regcoil_problem *handle);
int regcoil_c_setup(regcoil_problem *handle, const char *path);
int regcoil_c_build_matrices(regcoil_problem *handle);
int regcoil_c_prepare_solve(regcoil_problem *handle);
int regcoil_c_solve_ilambda(regcoil_problem *handle, int ilambda, double *chi2_b,
                            double *chi2_k, double *max_bn, double *max_k);
int regcoil_c_solve_lambda(regcoil_problem *handle, double lam, double *chi2_b,
                           double *chi2_k, double *max_bn, double *max_k);

#ifdef __cplusplus
}
#endif

#endif /* REGCOIL_C_API_H */
