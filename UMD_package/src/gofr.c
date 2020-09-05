/*
 * Author: Laurent Gilquin
 * Date:   September 2, 2020
 */ 

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <Python.h>
#include <numpy/arrayobject.h>
#include <numpy/npy_math.h>

#include "gofr.h"

static PyObject *compute_gofr_wrapper(PyObject *self, PyObject *args, 
                                      PyObject *kwargs) 
{
  PyArrayObject *X_, *res_, *coeff_, *types_;
  int mx, nx, nr, ntypes;
  double maxlength, discrete;
  long int *res;
  const long int *types;
  const double *X, *coeff;
  static char *kwlist[] = {"X", "res", "coeff", "types", "maxlength", "discrete",
  "ntypes", NULL};
  if (!PyArg_ParseTupleAndKeywords(args, kwargs, 
            "O!O!O!O!ddi:compute_gofr_wrapper", kwlist,
            &PyArray_Type, &X_,
            &PyArray_Type, &res_,
            &PyArray_Type, &coeff_,
            &PyArray_Type, &types_,
            &maxlength, &discrete, &ntypes)) {
    return NULL;
  }
  else {
    NPY_BEGIN_ALLOW_THREADS;
    X = PyArray_DATA(X_);
    res = PyArray_DATA(res_);
    coeff = PyArray_DATA(coeff_);
    types = PyArray_DATA(types_);
    mx = PyArray_DIM(X_, 0);
    nx = PyArray_DIM(X_, 1);
    nr = PyArray_DIM(res_, 1);

    compute_gofr(X, res, coeff, types, maxlength, discrete, ntypes, mx, nx, nr);
    NPY_END_ALLOW_THREADS;
  }

  Py_RETURN_NONE;
}


static PyObject *print_gofr_wrapper(PyObject *self, PyObject *args, 
                                    PyObject *kwargs) 
{
  PyArrayObject *gofr_, *types_;
  char* filename = NULL;
  int ng, ntypes;
  double maxlength, discrete, normalization;
  const long int *gofr, *types;
  static char *kwlist[] = {"gofr", "types", "ntypes", "maxlength", "discrete",
  "normalization", "filename", NULL};
  if (!PyArg_ParseTupleAndKeywords(args, kwargs, 
            "O!O!iddds:print_gofr_wrapper", kwlist,
            &PyArray_Type, &gofr_,
            &PyArray_Type, &types_,
            &ntypes,
            &maxlength, &discrete, &normalization,
            &filename)) {
    return NULL;
  }
  else {
    NPY_BEGIN_ALLOW_THREADS;
    gofr = PyArray_DATA(gofr_);
    types = PyArray_DATA(types_);
    ng = PyArray_DIM(gofr_, 1);
    
    FILE *fp = fopen(filename, "a");
    print_gofr(gofr, types, ntypes, maxlength, discrete, normalization, ng, fp);
    fclose(fp);
    NPY_END_ALLOW_THREADS;
  }

  Py_RETURN_NONE;
}


static PyObject *umd_pdist_wrapper(PyObject *self, PyObject *args, 
                                   PyObject *kwargs) 
{
  PyArrayObject *X_, *res_, *coeff_;
  int m, n;
  double *res;
  const double *X, *coeff;
  static char *kwlist[] = {"X", "res", "coeff", NULL};
  if (!PyArg_ParseTupleAndKeywords(args, kwargs, 
            "O!O!O!:umd_pdist_wrapper", kwlist,
            &PyArray_Type, &X_,
            &PyArray_Type, &res_,
            &PyArray_Type, &coeff_)) {
    return NULL;
  }
  else {
    NPY_BEGIN_ALLOW_THREADS;
    X = PyArray_DATA(X_);
    res = PyArray_DATA(res_);
    coeff = PyArray_DATA(coeff_);
    m = PyArray_DIM(X_, 0);
    n = PyArray_DIM(X_, 1);

    umd_distance(X, coeff, res, m, n);
    NPY_END_ALLOW_THREADS;
  }

  Py_RETURN_NONE;
}


// Method definition object for this extension, these argumens mean:
// ml_name: The name of the method
// ml_meth: Function pointer to the method implementation
// ml_flags: Flags indicating special features of this method, such as
//          accepting arguments, accepting keyword arguments, being a
//          class method, or being a static method of a class.
// ml_doc:  Contents of this method's docstring
static PyMethodDef c_gofr_method[] = {

    {"compute_gofr_wrapper",
    (PyCFunction) compute_gofr_wrapper,
    METH_VARARGS | METH_KEYWORDS,
    "To be filled."},
    {"print_gofr_wrapper",
    (PyCFunction) print_gofr_wrapper,
    METH_VARARGS | METH_KEYWORDS,
    "To be filled."},
    {"umd_pdist_wrapper",
    (PyCFunction) umd_pdist_wrapper,
    METH_VARARGS | METH_KEYWORDS,
    "To be filled."},
    {NULL, NULL, 0, NULL}        /* Sentinel */
};


// Module definition
// The arguments of this structure tell Python what to call your extension,
// what it's methods are and where to look for it's method definitions
static struct PyModuleDef c_gofr_definition = { 
    PyModuleDef_HEAD_INIT,
    "c_gofr",
    "To be filled.",
    -1, 
    c_gofr_method
};


// Module initialization
// Python calls this function when importing your extension. It is important
// that this function is named PyInit_[[your_module_name]] exactly, and matches
// the name keyword argument in setup.py's setup() call.
PyMODINIT_FUNC PyInit_c_gofr(void) {

    PyObject *m;
    
    m = PyModule_Create(&c_gofr_definition);
    import_array();
    
    return m;
}
