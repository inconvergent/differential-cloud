# -*- coding: utf-8 -*-
# cython: profile=True

from __future__ import division

cimport cloud
cimport numpy as np

from libc.stdlib cimport malloc, free
from zonemap3d cimport Zonemap3d


cdef class DifferentialCloud(cloud.Cloud):

  cdef double nearl
  cdef double midl
  cdef double farl

  cdef double *DX
  cdef double *DY
  cdef double *DZ

  cdef double *interactions
  cdef double *distances
  cdef long ptnum

  ## FUNCTIONS

  cdef double __get_interaction_by_types(self, long i, long j) nogil

  cdef double __get_interaction_by_id(self, long i, long j) nogil

  cdef double __get_particle_rad(self, long i) nogil

  cdef double __get_particles_radii(self, long i, long j) nogil

  cdef double __get_particles_harmonic_radii(self, long i, long j) nogil

  cdef double __get_particle_far(self, long i) nogil

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
    np.ndarray[double, mode="c",ndim=2] interactions,
    np.ndarray[double, mode="c",ndim=2] distances
  )

  cpdef long spawn(
    self,
    np.ndarray[long, mode="c",ndim=1] rnd
  )

  cpdef long optimize_position(
    self,
    double reject_stp,
    double attract_stp,
    double limit_stp
  )

