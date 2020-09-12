#include <stdio.h>
#include <Python.h>
#include <numpy/arrayobject.h>

double c_thetaneurons(){
  PyObject *obj;
  PyArrayIterObject *iter;
  iter = (PyArrayIterObject *)PyArray_IterNew(obj);
  while (iter->index < iter->size) {
    /* do something with the data at it->dataptr */
    PyArray_ITER_NEXT(iter);
  }
  double returnme = 1.0;
  return returnme;
}


// double c_thetaneurons(float t, x, e, KdivN, a):
//     # Model a network of oscillators
//     I_sync = cython_pulse(x).sum()
//     return (1 - cos(x)) + (1 + cos(x)) * (e + a * KdivN * I_sync);
