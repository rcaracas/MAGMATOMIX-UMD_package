/*
 * Author: Laurent Gilquin
 * Date:   August 26, 2020
 */ 

#if !defined(__clang__) && defined(__GNUC__) && defined(__GNUC_MINOR__)
#if __GNUC__ >= 5 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 4)
/* enable auto-vectorizer */
#pragma GCC optimize("tree-vectorize")
/* float associativity required to vectorize reductions */
#pragma GCC optimize("unsafe-math-optimizations")
/* maybe 5% gain, manual unrolling with more accumulators would be better */
#pragma GCC optimize("unroll-loops")
#endif
#endif

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <Python.h>
#include <numpy/arrayobject.h>
#include <numpy/npy_math.h>

#include "umd_dist.h"

static PyObject *umd_distance_wrapper(PyObject *self, PyObject *args, 
                                      PyObject *kwargs) 
{
  PyArrayObject *X_, *res_, *coeff_;
  int m, n;
  double *res;
  const double *X, *coeff;
  static char *kwlist[] = {"X", "res", "coeff", NULL};
  if (!PyArg_ParseTupleAndKeywords(args, kwargs, 
            "O!O!O!:umd_distance_wrapper", kwlist,
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
static PyMethodDef umd_distance_method[] = {

    {"umd_distance_wrapper",
    (PyCFunction) umd_distance_wrapper,
    METH_VARARGS | METH_KEYWORDS,
    "Calculate the umd distance, a shifted euclidean distance, in a C extension."},
    {NULL, NULL, 0, NULL}        /* Sentinel */
};



// Module definition
// The arguments of this structure tell Python what to call your extension,
// what it's methods are and where to look for it's method definitions
static struct PyModuleDef umd_distance_definition = { 
    PyModuleDef_HEAD_INIT,
    "umd_distance_wrapper",
    "A Python module that calculates the umd distance from C code.",
    -1, 
    umd_distance_method
};


// Module initialization
// Python calls this function when importing your extension. It is important
// that this function is named PyInit_[[your_module_name]] exactly, and matches
// the name keyword argument in setup.py's setup() call.
PyMODINIT_FUNC PyInit_umd_distance(void) {

    PyObject *m;
    
    m = PyModule_Create(&umd_distance_definition);
    import_array();
    
    return m;
}
