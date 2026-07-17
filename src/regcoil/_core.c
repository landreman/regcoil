#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include "regcoil_c_api.h"

static PyObject *
core_set_verbose(PyObject *self, PyObject *args)
{
    int flag;
    (void)self;
    if (!PyArg_ParseTuple(args, "p", &flag)) {
        return NULL;
    }
    regcoil_c_set_verbose(flag);
    Py_RETURN_NONE;
}

static PyObject *
core_nlambda(PyObject *self, PyObject *Py_UNUSED(ignored))
{
    (void)self;
    return PyLong_FromLong(regcoil_c_nlambda());
}

static PyObject *
core_setup(PyObject *self, PyObject *args)
{
    const char *path;
    int ierr;
    (void)self;
    if (!PyArg_ParseTuple(args, "s", &path)) {
        return NULL;
    }
    ierr = regcoil_c_setup(path);
    if (ierr != 0) {
        return PyErr_Format(PyExc_RuntimeError,
                            "regcoil setup failed with code %d for path %s",
                            ierr, path);
    }
    Py_RETURN_NONE;
}

static PyObject *
core_build_matrices(PyObject *self, PyObject *Py_UNUSED(ignored))
{
    int ierr;
    (void)self;
    ierr = regcoil_c_build_matrices();
    if (ierr != 0) {
        return PyErr_Format(PyExc_RuntimeError,
                            "regcoil_build_matrices failed with code %d", ierr);
    }
    Py_RETURN_NONE;
}

static PyObject *
core_prepare_solve(PyObject *self, PyObject *Py_UNUSED(ignored))
{
    int ierr;
    (void)self;
    ierr = regcoil_c_prepare_solve();
    if (ierr != 0) {
        return PyErr_Format(PyExc_RuntimeError,
                            "regcoil_prepare_solve failed with code %d", ierr);
    }
    Py_RETURN_NONE;
}

static PyObject *
metrics_tuple(double chi2_b, double chi2_k, double max_bn, double max_k)
{
    return Py_BuildValue("dddd", chi2_b, chi2_k, max_bn, max_k);
}

static PyObject *
core_solve_ilambda(PyObject *self, PyObject *args)
{
    int ilambda0;
    int ierr;
    double chi2_b, chi2_k, max_bn, max_k;
    (void)self;
    if (!PyArg_ParseTuple(args, "i", &ilambda0)) {
        return NULL;
    }
    /* Python uses 0-based indices; Fortran uses 1-based. */
    ierr = regcoil_c_solve_ilambda(ilambda0 + 1, &chi2_b, &chi2_k, &max_bn, &max_k);
    if (ierr != 0) {
        return PyErr_Format(PyExc_RuntimeError,
                            "regcoil solve_ilambda(%d) failed with code %d",
                            ilambda0, ierr);
    }
    return metrics_tuple(chi2_b, chi2_k, max_bn, max_k);
}

static PyObject *
core_solve_lambda(PyObject *self, PyObject *args)
{
    double lam;
    int ierr;
    double chi2_b, chi2_k, max_bn, max_k;
    (void)self;
    if (!PyArg_ParseTuple(args, "d", &lam)) {
        return NULL;
    }
    ierr = regcoil_c_solve_lambda(lam, &chi2_b, &chi2_k, &max_bn, &max_k);
    if (ierr != 0) {
        return PyErr_Format(PyExc_RuntimeError,
                            "regcoil solve_lambda(%g) failed with code %d",
                            lam, ierr);
    }
    return metrics_tuple(chi2_b, chi2_k, max_bn, max_k);
}

static PyMethodDef core_methods[] = {
    {"set_verbose", core_set_verbose, METH_VARARGS,
     "set_verbose(flag)\n\nEnable or disable Fortran verbose logging."},
    {"nlambda", core_nlambda, METH_NOARGS,
     "nlambda()\n\nReturn the number of lambda values from the last setup()."},
    {"setup", core_setup, METH_VARARGS,
     "setup(path)\n\nRead namelist, build matrices, and prepare the solve (globals)."},
    {"build_matrices", core_build_matrices, METH_NOARGS,
     "build_matrices()\n\nAssemble REGCOIL matrices (requires prior geometry init)."},
    {"prepare_solve", core_prepare_solve, METH_NOARGS,
     "prepare_solve()\n\nAllocate work arrays for regularized solves."},
    {"solve_ilambda", core_solve_ilambda, METH_VARARGS,
     "solve_ilambda(i)\n\nSolve at lambda index i (0-based). Returns "
     "(chi2_B, chi2_K, max_Bnormal, max_K)."},
    {"solve_lambda", core_solve_lambda, METH_VARARGS,
     "solve_lambda(lambda_)\n\nSolve at the given lambda (uses slot 0). Returns "
     "(chi2_B, chi2_K, max_Bnormal, max_K)."},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef core_module = {
    PyModuleDef_HEAD_INIT,
    .m_name = "regcoil._core",
    .m_doc = "Fortran REGCOIL core extension (iso_c_binding).",
    .m_size = -1,
    .m_methods = core_methods,
};

PyMODINIT_FUNC
PyInit__core(void)
{
    return PyModule_Create(&core_module);
}
