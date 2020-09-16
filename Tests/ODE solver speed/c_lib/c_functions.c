#include <stdio.h>
#include <math.h>
#include <Python.h>
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>

// Tried an approach like https://github.com/johnnylee/python-numpy-c-extension-examples
// Functions:
static PyObject *try_print(PyObject* self, PyObject *args){
  int pid, sts=0;
  if(!PyArg_ParseTuple(args, "i", &pid))
  {
    return NULL;
  }
  printf("Hello, from C World! Pid: %i \n", pid);
  sts=pid;
  //return Py_BuildValue("i", sts);
  Py_RETURN_NONE;
}

double c_factorial(int n) {
    if (n==1 || n==0) return 1;
    return n * c_factorial(n-1);
}

static PyObject *factorial(PyObject *self, PyObject *args){
  int m = 100, *n = &m;
  printf("Hello, from C World! Pid: %i \n", *n);

  if(!PyArg_ParseTuple(args, "i", &n)) {
        return NULL;
    }
  double result = c_factorial(*n);
  return Py_BuildValue("i", result);
}

// Python structs to make it compile
static PyMethodDef c_methods[] = {
    {"printing", try_print, METH_VARARGS, "print function."},
    {"factorial", factorial, METH_VARARGS, "factorial function."},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef c_functions_module = {
    PyModuleDef_HEAD_INIT, "c_functions", "Extra C functions.", -1, c_methods
};

PyMODINIT_FUNC PyInit_c_functions(void) {
    import_array();
    return PyModule_Create(&c_functions_module);
}
