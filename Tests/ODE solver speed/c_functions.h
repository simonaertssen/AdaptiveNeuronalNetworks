// Using the numpy c API see
// https://stackoverflow.com/questions/56182259/how-does-one-acces-numpy-multidimensionnal-array-in-c-extensions
// https://scipy-lectures.org/advanced/interfacing_with_c/interfacing_with_c.html

#include <stdio.h>
#include <math.h>
#include <Python.h>
#include <numpy/arrayobject.h>

// Pure C functions
double factorial(int n){
  if (n==1 || n==0) return 1;
  else return n * factorial(n-1);
}

double c_a_n(int n){
  return pow(2,n)*pow(factorial(n),2)/factorial(2*n);
}


// Forward numpy declarations:
static PyObject *c_pulse(PyObject *self, PyObject *args);
static PyMethodDef Cmethods[] = {
     {"c_pulse", c_pulse, METH_VARARGS, "Make pulse-shaped cosine wave."},
     {NULL, NULL, 0, NULL}
};

static struct PyModuleDef cModPyDem = {
    PyModuleDef_HEAD_INIT,
    "c_functions", "Speeding up a DOPRI54 ODE solver.", -1, Cmethods
};

PyMODINIT_FUNC PyInit_cos_module(void) {
    return PyModule_Create(&cModPyDem);
}


// Function declarations:
static PyObject *c_pulse(PyObject *theta):
    return theta

// static PyObject *c_thetaneurons(){
//
//   I_sync = numpy_pulse(x).sum()
//   return (1 - cos(x)) + (1 + cos(x)) * (e + a * KdivN * I_sync);
// }
