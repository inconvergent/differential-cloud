# -*- coding: utf-8 -*-
# cython: profile=True

from __future__ import division

cimport cloud
cimport numpy as np

from libc.stdlib cimport malloc, free
from zonemap3d cimport Zonemap3d


cdef class DifferentialCloud(cloud.Cloud):

  cdef double *DX
  cdef double *DY
  cdef double *DZ

  cdef long num_particle_types
  cdef double *interactions
  cdef double *distances

  ## FUNCTIONS

  cdef double __get_interaction_by_types(self, long i, long j) nogil

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

  cdef double __get_particle_rad(
    self,
    long
  ) nogil

  cdef double __get_particle_far(
    self,
    long
  ) nogil

  cdef double __get_particles_radii(self, long i, long j) nogil

  cdef double __get_particles_harmonic_radii(self, long i, long j) nogil

  ## EXTERNAL

  cpdef long position_noise(
    self,
    np.ndarray[double, mode="c",ndim=2] a
  )

  cpdef long init_rules(
    self,
    long num_particle_types,
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
    double stp_limit
  )

