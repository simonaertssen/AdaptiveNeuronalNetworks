#include <stdio.h>
#include <math.h>
#include <Python.h>
#include <numpy/arrayobject.h>

static PyObject *intsumtest(PyObject *self, PyObject *args){
  int *a, *b;
  if (!PyArg_ParseTuple(args, "ii", &a, &b)) {
    return NULL;
  }
  long int sum = (*a) + (*b);
  printf("%ld + %ld = %ld", *a, *b, sum);
  //return Py_BuildValue("i", sum);
}


static PyObject *sumtest(PyObject *self, PyObject *args){
  PyArrayObject *np_theta;
  npy_float64 *theta;
  npy_int64 length;
  if (!PyArg_ParseTuple(args, "lO!", &length, &PyArray_Type, &np_theta)) {
    return NULL;
  }

  theta = (npy_float64*)PyArray_DATA(np_theta);

  long i = 0;
  double sum = 0;
  for (i = 0; i < length; ++i){
    sum += theta[i];
  }
  return PyFloat_FromDouble(sum);
}
