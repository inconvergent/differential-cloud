# -*- coding: utf-8 -*-
# cython: profile=True

from __future__ import division

cimport cloud
cimport numpy as np

from libc.stdlib cimport malloc, free
from zonemap3d cimport Zonemap3d


cdef class DifferentialCloud(cloud.Cloud):

  cdef double nearl
  cdef double farl

  cdef double *DX
  cdef double *DY
  cdef double *DZ

  cdef long *rules
  cdef long num_rules

  ## FUNCTIONS

  cdef long __get_rule_by_types(self, long i, long j) nogil

  cdef long __rules(
    self,
    long v1,
    double *dx,
    double *dy,
    double *dz,
    double stp,
    long *vertices,
    double *dst,
    long num
  ) nogil

  ## EXTERNAL

  cpdef long position_noise(
    self,
    np.ndarray[double, mode="c",ndim=2] a
  )

  cpdef long init_rules(
    self,
    np.ndarray[long, mode="c",ndim=2] r
  )

  cpdef long optimize_position(
    self,
    double reject_stp,
    double attract_stp
  )

