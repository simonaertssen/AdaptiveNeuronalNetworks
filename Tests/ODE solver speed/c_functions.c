#include <stdio.h>
#include <math.h>
#include <Python.h>
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>


// Forward declarations:
static PyObject *factorial(PyObject *n, PyObject *args);
long c_factorial(long *n);

static PyObject *c_a_n(PyObject *n, PyObject *args);

// Python structs to make it compile
static PyMethodDef c_methods[] = {
    {"factorial", factorial, METH_VARARGS, "factorial function."},
    {"c_a_n", c_a_n, METH_VARARGS, "a_n function."},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef c_functions_module = {
    PyModuleDef_HEAD_INIT, "c_functions", "Extra C functions.", -1, c_methods
};

PyMODINIT_FUNC PyInit_c_functions(void) {
    import_array();
    return PyModule_Create(&c_functions_module);
}


// Functions:
static PyObject *factorial(PyObject *self, PyObject *args){
  long *n;
  if (!PyArg_ParseTuple(args, "l", &n)) {
        return NULL;
    }
  return Py_BuildValue("l", c_factorial(n));
}

long c_factorial(long *n) {
    if (*n==1 || *n==0) return 1;
    return (long)n * c_factorial(n-1);
}

static PyObject *c_a_n(PyObject *self, PyObject *args){
  long *n, *m;
  if (!PyArg_ParseTuple(args, "l", &n)) {
        return NULL;
  }
  m = &(2*(*n));
  return Py_BuildValue("l", pow(2,(double)(*n))*pow(c_factorial(n),2)/c_factorial(m));
}
