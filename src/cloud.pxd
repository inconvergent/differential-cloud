# -*- coding: utf-8 -*-
# cython: profile=True

cimport numpy as np
from zonemap3d cimport Zonemap3d

cdef class Cloud:

  cdef long nmax
  cdef long vnum

  cdef long nz
  cdef double zonewidth

  cdef long procs
  cdef long state
  cdef double start_time

  ## ARRAYS

  cdef double *X
  cdef double *Y
  cdef double *Z
  cdef long *A

  cdef double *rules

  cdef Zonemap3d zonemap

  ## FUNCTIONS

  ## INTERNAL

  cdef long __valid_new_vertex(self, double x, double y, double z) nogil

  cdef long __new_vertex(self, double x, double y, double z) nogil

  cdef void __set_vertex_type(self, long v1, long i) nogil

  ## EXTERNAL

  cpdef long np_get_vertices(self, np.ndarray[double, mode="c",ndim=2] x)

  cpdef long init_cloud(
    self,
    np.ndarray[double, mode="c",ndim=2] x,
    np.ndarray[long, mode="c",ndim=1] a
  )

  cpdef long set_rules(
    self,
    np.ndarray[double, mode="c",ndim=2] r
  )

  ## INFO

  cpdef long get_vnum(self)

  cpdef double get_start_time(self)

