#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include "regcoil_c_api.h"

typedef struct {
    PyObject_HEAD
    regcoil_problem *handle;
} RegcoilProblemObject;

static int
RegcoilProblem_clear_handle(RegcoilProblemObject *self)
{
    if (self->handle != NULL) {
        regcoil_c_destroy(self->handle);
        self->handle = NULL;
    }
    return 0;
}

static void
RegcoilProblem_dealloc(RegcoilProblemObject *self)
{
    RegcoilProblem_clear_handle(self);
    Py_TYPE(self)->tp_free((PyObject *)self);
}

static int
RegcoilProblem_init(RegcoilProblemObject *self, PyObject *args, PyObject *kw)
{
    (void)args;
    (void)kw;
    if (self->handle != NULL) {
        RegcoilProblem_clear_handle(self);
    }
    self->handle = regcoil_c_create();
    if (self->handle == NULL) {
        PyErr_SetString(PyExc_MemoryError, "regcoil_c_create failed");
        return -1;
    }
    return 0;
}

static PyObject *
RegcoilProblem_set_verbose(RegcoilProblemObject *self, PyObject *args)
{
    int flag;
    if (!PyArg_ParseTuple(args, "p", &flag)) {
        return NULL;
    }
    regcoil_c_set_verbose(self->handle, flag);
    Py_RETURN_NONE;
}

static PyObject *
RegcoilProblem_nlambda(RegcoilProblemObject *self, PyObject *Py_UNUSED(ignored))
{
    return PyLong_FromLong(regcoil_c_nlambda(self->handle));
}

static PyObject *
RegcoilProblem_setup(RegcoilProblemObject *self, PyObject *args)
{
    const char *path;
    int ierr;
    if (!PyArg_ParseTuple(args, "s", &path)) {
        return NULL;
    }
    ierr = regcoil_c_setup(self->handle, path);
    if (ierr != 0) {
        return PyErr_Format(PyExc_RuntimeError,
                            "regcoil setup failed with code %d for path %s",
                            ierr, path);
    }
    Py_RETURN_NONE;
}

static PyObject *
RegcoilProblem_build_matrices(RegcoilProblemObject *self,
                              PyObject *Py_UNUSED(ignored))
{
    int ierr = regcoil_c_build_matrices(self->handle);
    if (ierr != 0) {
        return PyErr_Format(PyExc_RuntimeError,
                            "regcoil_build_matrices failed with code %d", ierr);
    }
    Py_RETURN_NONE;
}

static PyObject *
RegcoilProblem_prepare_solve(RegcoilProblemObject *self,
                             PyObject *Py_UNUSED(ignored))
{
    int ierr = regcoil_c_prepare_solve(self->handle);
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
RegcoilProblem_solve_ilambda(RegcoilProblemObject *self, PyObject *args)
{
    int ilambda0;
    int ierr;
    double chi2_b, chi2_k, max_bn, max_k;
    if (!PyArg_ParseTuple(args, "i", &ilambda0)) {
        return NULL;
    }
    ierr = regcoil_c_solve_ilambda(self->handle, ilambda0 + 1, &chi2_b, &chi2_k,
                                   &max_bn, &max_k);
    if (ierr != 0) {
        return PyErr_Format(PyExc_RuntimeError,
                            "regcoil solve_ilambda(%d) failed with code %d",
                            ilambda0, ierr);
    }
    return metrics_tuple(chi2_b, chi2_k, max_bn, max_k);
}

static PyObject *
RegcoilProblem_solve_lambda(RegcoilProblemObject *self, PyObject *args)
{
    double lam;
    int ierr;
    double chi2_b, chi2_k, max_bn, max_k;
    if (!PyArg_ParseTuple(args, "d", &lam)) {
        return NULL;
    }
    ierr = regcoil_c_solve_lambda(self->handle, lam, &chi2_b, &chi2_k, &max_bn,
                                  &max_k);
    if (ierr != 0) {
        return PyErr_Format(PyExc_RuntimeError,
                            "regcoil solve_lambda(%g) failed with code %d", lam,
                            ierr);
    }
    return metrics_tuple(chi2_b, chi2_k, max_bn, max_k);
}

static PyMethodDef RegcoilProblem_methods[] = {
    {"set_verbose", (PyCFunction)RegcoilProblem_set_verbose, METH_VARARGS,
     "set_verbose(flag)"},
    {"nlambda", (PyCFunction)RegcoilProblem_nlambda, METH_NOARGS, "nlambda()"},
    {"setup", (PyCFunction)RegcoilProblem_setup, METH_VARARGS,
     "setup(path)\n\nRead namelist, build matrices, prepare solve."},
    {"build_matrices", (PyCFunction)RegcoilProblem_build_matrices, METH_NOARGS,
     "build_matrices()"},
    {"prepare_solve", (PyCFunction)RegcoilProblem_prepare_solve, METH_NOARGS,
     "prepare_solve()"},
    {"solve_ilambda", (PyCFunction)RegcoilProblem_solve_ilambda, METH_VARARGS,
     "solve_ilambda(i) -> (chi2_B, chi2_K, max_Bnormal, max_K)"},
    {"solve_lambda", (PyCFunction)RegcoilProblem_solve_lambda, METH_VARARGS,
     "solve_lambda(lambda_) -> (chi2_B, chi2_K, max_Bnormal, max_K)"},
    {NULL, NULL, 0, NULL},
};

static PyTypeObject RegcoilProblemType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    .tp_name = "regcoil._core.RegcoilProblem",
    .tp_basicsize = sizeof(RegcoilProblemObject),
    .tp_dealloc = (destructor)RegcoilProblem_dealloc,
    .tp_flags = Py_TPFLAGS_DEFAULT,
    .tp_doc = "Opaque handle to one Fortran REGCOIL problem instance.",
    .tp_methods = RegcoilProblem_methods,
    .tp_init = (initproc)RegcoilProblem_init,
    .tp_new = PyType_GenericNew,
};

static struct PyModuleDef core_module = {
    PyModuleDef_HEAD_INIT,
    .m_name = "regcoil._core",
    .m_doc = "Fortran REGCOIL core extension (instance handles).",
    .m_size = -1,
};

PyMODINIT_FUNC
PyInit__core(void)
{
    PyObject *m;
    if (PyType_Ready(&RegcoilProblemType) < 0) {
        return NULL;
    }
    m = PyModule_Create(&core_module);
    if (m == NULL) {
        return NULL;
    }
    Py_INCREF(&RegcoilProblemType);
    if (PyModule_AddObject(m, "RegcoilProblem", (PyObject *)&RegcoilProblemType) <
        0) {
        Py_DECREF(&RegcoilProblemType);
        Py_DECREF(m);
        return NULL;
    }
    return m;
}
