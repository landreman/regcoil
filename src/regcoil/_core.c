#define PY_SSIZE_T_CLEAN
#include <Python.h>
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>
#include "regcoil_c_api.h"

/* Module-level functions wrapping the stateless Fortran kernels.
 * Large arrays (r_plasma, r_coil, basis_functions, ...) must already be
 * float64 and Fortran-contiguous: a silent copy of a multi-MB array would
 * be a hidden cost, so we raise an error instead of casting (see API.md conventions). 
 * */

static PyArrayObject *
require_f64_farray(PyObject *obj, int ndim, const char *name)
{
    PyArrayObject *arr;
    if (!PyArray_Check(obj)) {
        PyErr_Format(PyExc_TypeError, "%s must be a numpy array", name);
        return NULL;
    }
    arr = (PyArrayObject *)obj;
    if (PyArray_TYPE(arr) != NPY_DOUBLE) {
        PyErr_Format(PyExc_TypeError, "%s must have dtype float64", name);
        return NULL;
    }
    if (PyArray_NDIM(arr) != ndim) {
        PyErr_Format(PyExc_ValueError, "%s must have ndim=%d, got %d", name, ndim, PyArray_NDIM(arr));
        return NULL;
    }
    if (ndim > 1 && !PyArray_IS_F_CONTIGUOUS(arr)) {
        PyErr_Format(PyExc_ValueError,
                     "%s must be Fortran-contiguous (pass order='F' or np.asfortranarray(...))", name);
        return NULL;
    }
    return arr;
}

static int
dims_match(PyArrayObject *a, PyArrayObject *b)
{
    int i;
    if (PyArray_NDIM(a) != PyArray_NDIM(b)) {
        return 0;
    }
    for (i = 0; i < PyArray_NDIM(a); i++) {
        if (PyArray_DIM(a, i) != PyArray_DIM(b, i)) {
            return 0;
        }
    }
    return 1;
}

static PyObject *
core_omp_max_threads(PyObject *Py_UNUSED(self), PyObject *Py_UNUSED(ignored))
{
    return PyLong_FromLong((long)regcoil_c_omp_max_threads());
}

static PyObject *
core_build_inductance(PyObject *Py_UNUSED(self), PyObject *args)
{
    PyObject *r_plasma_o, *normal_plasma_o, *r_coil_o, *normal_coil_o, *drdtheta_coil_o, *drdzeta_coil_o;
    int nfp;
    double net_poloidal_current, net_toroidal_current, dtheta_coil, dzeta_coil;
    PyArrayObject *r_plasma, *normal_plasma, *r_coil, *normal_coil, *drdtheta_coil, *drdzeta_coil;
    npy_intp nt_p, nz_p, nt_c, nzetal_c, nz_c;
    npy_intp inductance_dims[2], h_dims[1];
    PyArrayObject *inductance, *h;
    int ierr;

    if (!PyArg_ParseTuple(args, "OOOOOOidddd", &r_plasma_o, &normal_plasma_o, &r_coil_o, &normal_coil_o,
                          &drdtheta_coil_o, &drdzeta_coil_o, &nfp, &net_poloidal_current,
                          &net_toroidal_current, &dtheta_coil, &dzeta_coil)) {
        return NULL;
    }
    if (nfp < 1) {
        PyErr_SetString(PyExc_ValueError, "nfp must be a positive integer");
        return NULL;
    }

    r_plasma = require_f64_farray(r_plasma_o, 3, "r_plasma");
    normal_plasma = require_f64_farray(normal_plasma_o, 3, "normal_plasma");
    r_coil = require_f64_farray(r_coil_o, 3, "r_coil");
    normal_coil = require_f64_farray(normal_coil_o, 3, "normal_coil");
    drdtheta_coil = require_f64_farray(drdtheta_coil_o, 3, "drdtheta_coil");
    drdzeta_coil = require_f64_farray(drdzeta_coil_o, 3, "drdzeta_coil");
    if (!r_plasma || !normal_plasma || !r_coil || !normal_coil || !drdtheta_coil || !drdzeta_coil) {
        return NULL;
    }
    if (!dims_match(r_plasma, normal_plasma)) {
        PyErr_SetString(PyExc_ValueError, "r_plasma and normal_plasma must have the same shape");
        return NULL;
    }
    if (!dims_match(r_coil, normal_coil) || !dims_match(r_coil, drdtheta_coil) || !dims_match(r_coil, drdzeta_coil)) {
        PyErr_SetString(PyExc_ValueError, "r_coil, normal_coil, drdtheta_coil, drdzeta_coil must have the same shape");
        return NULL;
    }
    nt_p = PyArray_DIM(r_plasma, 1);
    nz_p = PyArray_DIM(r_plasma, 2);
    nt_c = PyArray_DIM(r_coil, 1);
    nzetal_c = PyArray_DIM(r_coil, 2);
    if (nzetal_c % nfp != 0) {
        PyErr_SetString(PyExc_ValueError, "r_coil's last extent must be a multiple of nfp");
        return NULL;
    }
    nz_c = nzetal_c / nfp;

    inductance_dims[0] = nt_p * nz_p;
    inductance_dims[1] = nt_c * nz_c;
    inductance = (PyArrayObject *)PyArray_ZEROS(2, inductance_dims, NPY_DOUBLE, 1 /* fortran order */);
    if (!inductance) {
        return NULL;
    }
    h_dims[0] = nt_p * nz_p;
    h = (PyArrayObject *)PyArray_ZEROS(1, h_dims, NPY_DOUBLE, 0);
    if (!h) {
        Py_DECREF(inductance);
        return NULL;
    }

    Py_BEGIN_ALLOW_THREADS
    ierr = regcoil_c_build_inductance(
        (int)nt_p, (int)nz_p, (int)nt_c, (int)nz_c, nfp,
        (const double *)PyArray_DATA(r_plasma), (const double *)PyArray_DATA(normal_plasma),
        (const double *)PyArray_DATA(r_coil), (const double *)PyArray_DATA(normal_coil),
        (const double *)PyArray_DATA(drdtheta_coil), (const double *)PyArray_DATA(drdzeta_coil),
        net_poloidal_current, net_toroidal_current, dtheta_coil, dzeta_coil,
        (double *)PyArray_DATA(inductance), (double *)PyArray_DATA(h));
    Py_END_ALLOW_THREADS
    if (ierr != 0) {
        Py_DECREF(inductance);
        Py_DECREF(h);
        return PyErr_Format(PyExc_RuntimeError, "regcoil_build_inductance failed with info=%d", ierr);
    }

    return Py_BuildValue("NN", inductance, h);
}

static PyObject *
core_build_g_and_h(PyObject *Py_UNUSED(self), PyObject *args)
{
    PyObject *r_plasma_o, *normal_plasma_o, *r_coil_o, *normal_coil_o, *drdtheta_coil_o, *drdzeta_coil_o;
    PyObject *basis_functions_o;
    int nfp;
    double net_poloidal_current, net_toroidal_current, dtheta_coil, dzeta_coil;
    PyArrayObject *r_plasma, *normal_plasma, *r_coil, *normal_coil, *drdtheta_coil, *drdzeta_coil;
    PyArrayObject *basis_functions;
    npy_intp nt_p, nz_p, nt_c, nzetal_c, nz_c, nbf;
    npy_intp g_dims[2], h_dims[1];
    PyArrayObject *g, *h;
    int ierr;

    if (!PyArg_ParseTuple(args, "OOOOOOOidddd", &r_plasma_o, &normal_plasma_o, &r_coil_o, &normal_coil_o,
                          &drdtheta_coil_o, &drdzeta_coil_o, &basis_functions_o, &nfp,
                          &net_poloidal_current, &net_toroidal_current, &dtheta_coil, &dzeta_coil)) {
        return NULL;
    }
    if (nfp < 1) {
        PyErr_SetString(PyExc_ValueError, "nfp must be a positive integer");
        return NULL;
    }

    r_plasma = require_f64_farray(r_plasma_o, 3, "r_plasma");
    normal_plasma = require_f64_farray(normal_plasma_o, 3, "normal_plasma");
    r_coil = require_f64_farray(r_coil_o, 3, "r_coil");
    normal_coil = require_f64_farray(normal_coil_o, 3, "normal_coil");
    drdtheta_coil = require_f64_farray(drdtheta_coil_o, 3, "drdtheta_coil");
    drdzeta_coil = require_f64_farray(drdzeta_coil_o, 3, "drdzeta_coil");
    basis_functions = require_f64_farray(basis_functions_o, 2, "basis_functions");
    if (!r_plasma || !normal_plasma || !r_coil || !normal_coil || !drdtheta_coil || !drdzeta_coil || !basis_functions) {
        return NULL;
    }
    if (!dims_match(r_plasma, normal_plasma)) {
        PyErr_SetString(PyExc_ValueError, "r_plasma and normal_plasma must have the same shape");
        return NULL;
    }
    if (!dims_match(r_coil, normal_coil) || !dims_match(r_coil, drdtheta_coil) || !dims_match(r_coil, drdzeta_coil)) {
        PyErr_SetString(PyExc_ValueError, "r_coil, normal_coil, drdtheta_coil, drdzeta_coil must have the same shape");
        return NULL;
    }
    nt_p = PyArray_DIM(r_plasma, 1);
    nz_p = PyArray_DIM(r_plasma, 2);
    nt_c = PyArray_DIM(r_coil, 1);
    nzetal_c = PyArray_DIM(r_coil, 2);
    if (nzetal_c % nfp != 0) {
        PyErr_SetString(PyExc_ValueError, "r_coil's last extent must be a multiple of nfp");
        return NULL;
    }
    nz_c = nzetal_c / nfp;
    nbf = PyArray_DIM(basis_functions, 1);
    if (PyArray_DIM(basis_functions, 0) != nt_c * nz_c) {
        PyErr_SetString(PyExc_ValueError, "basis_functions must have shape (ntheta_coil*nzeta_coil, nbf)");
        return NULL;
    }

    g_dims[0] = nt_p * nz_p;
    g_dims[1] = nbf;
    g = (PyArrayObject *)PyArray_ZEROS(2, g_dims, NPY_DOUBLE, 1 /* fortran order */);
    if (!g) {
        return NULL;
    }
    h_dims[0] = nt_p * nz_p;
    h = (PyArrayObject *)PyArray_ZEROS(1, h_dims, NPY_DOUBLE, 0);
    if (!h) {
        Py_DECREF(g);
        return NULL;
    }

    Py_BEGIN_ALLOW_THREADS
    ierr = regcoil_c_build_g_and_h(
        (int)nt_p, (int)nz_p, (int)nt_c, (int)nz_c, nfp, (int)nbf,
        (const double *)PyArray_DATA(r_plasma), (const double *)PyArray_DATA(normal_plasma),
        (const double *)PyArray_DATA(r_coil), (const double *)PyArray_DATA(normal_coil),
        (const double *)PyArray_DATA(drdtheta_coil), (const double *)PyArray_DATA(drdzeta_coil),
        (const double *)PyArray_DATA(basis_functions),
        net_poloidal_current, net_toroidal_current, dtheta_coil, dzeta_coil,
        (double *)PyArray_DATA(g), (double *)PyArray_DATA(h));
    Py_END_ALLOW_THREADS
    if (ierr != 0) {
        Py_DECREF(g);
        Py_DECREF(h);
        return PyErr_Format(PyExc_RuntimeError, "regcoil_build_g_and_h failed with info=%d", ierr);
    }

    return Py_BuildValue("NN", g, h);
}

static PyMethodDef core_methods[] = {
    {"omp_max_threads", core_omp_max_threads, METH_NOARGS,
     "omp_max_threads() -> int\n\n"
     "Return omp_get_max_threads() from the compiled extension. This function is useful for determining how many threads are seen by the compiled kernels."},
    {"build_inductance", core_build_inductance, METH_VARARGS,
     "build_inductance(r_plasma, normal_plasma, r_coil, normal_coil, drdtheta_coil, drdzeta_coil,\n"
     "                  nfp, net_poloidal_current, net_toroidal_current, dtheta_coil, dzeta_coil)\n"
     "-> (inductance, h)\n\n"
     "Debug/regression variant: materializes the full inductance matrix."},
    {"build_g_and_h", core_build_g_and_h, METH_VARARGS,
     "build_g_and_h(r_plasma, normal_plasma, r_coil, normal_coil, drdtheta_coil, drdzeta_coil,\n"
     "              basis_functions, nfp, net_poloidal_current, net_toroidal_current,\n"
     "              dtheta_coil, dzeta_coil) -> (g, h)\n\n"
     "Fused kernel: g = dtheta_coil*dzeta_coil * (inductance @ basis_functions), without\n"
     "materializing the full inductance matrix."},
    {NULL, NULL, 0, NULL},
};

static struct PyModuleDef core_module = {
    PyModuleDef_HEAD_INIT,
    .m_name = "regcoil._core",
    .m_doc = "Fortran REGCOIL core extension",
    .m_size = -1,
    .m_methods = core_methods,
};

PyMODINIT_FUNC
PyInit__core(void)
{
    import_array();
    return PyModule_Create(&core_module);
}
