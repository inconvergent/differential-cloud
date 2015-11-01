# -*- coding: utf-8 -*-

cimport cython

cdef inline void long_array_init(long *a, long n, long v) nogil:
  cdef long i
  for i in xrange(n):
    a[i] = v
  return

cdef inline void int_array_init(int *a, long n, int v) nogil:
  cdef long i
  for i in xrange(n):
    a[i] = v
  return

cdef inline void float_array_init(float *a, long n, float v) nogil:
  cdef long i
  for i in xrange(n):
    a[i] = v
  return

cdef inline void double_array_init(double *a, long n, double v) nogil:
  cdef long i
  for i in xrange(n):
    a[i] = v
  return

