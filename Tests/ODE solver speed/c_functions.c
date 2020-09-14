#include <stdio.h>
#include <math.h>
#include <Python.h>
#include <numpy/arrayobject.h>

// Forward declarations:
double factorial(int n);
double c_a_n(int n);

// Python structs to make it compile
static PyMethodDef c_methods[] = {
    {"factorial", factorial, METH_VARARGS, "factorial function."},
    {"c_a_n", c_a_n, METH_VARARGS, "a_n function."},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef c_functions_module = {
    PyModuleDef_HEAD_INIT, "c_functions", "Extra C functions.", -1, c_methods
};

PyMODINIT_FUNC
PyInit_cos_module(void)
{
    import_array();
    return PyModule_Create(&c_functions_module);
}


// Functions:
double factorial(int n){
  if (n==1 || n==0) return 1;
  else return n * factorial(n-1);
}

double c_a_n(int n){
  return pow(2,n)*pow(factorial(n),2)/factorial(2*n);
}
