# -*- coding: utf-8 -*-
# cython: profile=True

from __future__ import division

cimport cython
cimport cloud

from zonemap3d cimport Zonemap3d

from cython.parallel import parallel, prange

from libc.math cimport sqrt
from libc.math cimport pow
from libc.math cimport cos
from libc.math cimport sin
from libc.math cimport fabs

from helpers cimport double_array_init
from helpers cimport long_array_init

import numpy as np
cimport numpy as np


cdef class DifferentialCloud(cloud.Cloud):

  def __init__(
    self,
    long nmax,
    double zonewidth,
    long procs
  ):

    cloud.Cloud.__init__(self, nmax, zonewidth, procs)

    self.num_particle_types = 0

    return

  def __cinit__(self, long nmax, *arg, **args):

    self.DX = <double *>malloc(nmax*sizeof(double))
    self.DY = <double *>malloc(nmax*sizeof(double))
    self.DZ = <double *>malloc(nmax*sizeof(double))

    return

  def __dealloc__(self):

    free(self.DX)
    free(self.DY)
    free(self.DZ)

    return

  @cython.wraparound(False)
  @cython.boundscheck(False)
  @cython.nonecheck(False)
  cdef double __get_interaction_by_types(self, long i, long j) nogil:

    return self.interactions[i*self.num_particle_types+j]

  @cython.wraparound(False)
  @cython.boundscheck(False)
  @cython.nonecheck(False)
  cdef double __get_particle_rad(
    self,
    long i
  ) nogil:

    return self.distances[i*2]

  @cython.wraparound(False)
  @cython.boundscheck(False)
  @cython.nonecheck(False)
  cdef double __get_particle_far(
    self,
    long i
  ) nogil:

    return self.distances[i*2+1]

  @cython.wraparound(False)
  @cython.boundscheck(False)
  @cython.nonecheck(False)
  @cython.cdivision(True)
  cdef double __get_particles_radii(self, long i, long j) nogil:

    return self.__get_particle_rad(i) + self.__get_particle_rad(j)

  @cython.wraparound(False)
  @cython.boundscheck(False)
  @cython.nonecheck(False)
  @cython.cdivision(True)
  cdef double __get_particles_harmonic_radii(self, long i, long j) nogil:

    cdef double ri = self.__get_particle_rad(i)
    cdef double rj = self.__get_particle_rad(j)

    return ri*rj/(ri+rj)

  @cython.wraparound(False)
  @cython.boundscheck(False)
  @cython.nonecheck(False)
  @cython.cdivision(True)
  cdef long __rules(
    self,
    long v,
    double *diffx,
    double *diffy,
    double *diffz,
    double stp,
    long *vertices,
    double *dst,
    long num
  ) nogil:

    cdef long k
    cdef long i
    cdef long k4

    cdef long tn
    cdef long tv

    cdef double dx
    cdef double dy
    cdef double dz
    cdef double nrm
    cdef double s

    cdef double rule
    cdef double radii

    cdef double resx = 0.0
    cdef double resy = 0.0
    cdef double resz = 0.0

    cdef double vrad = self.__get_particle_rad(v)
    cdef double vfar = self.__get_particle_far(v)

    for k in range(num):

      neigh = vertices[k]

      if neigh == v:
        continue

      k4 = k*4
      nrm = dst[k4+3]

      if nrm < 1.e-10 or nrm>vfar:
        continue

      tv = self.A[v]
      tn = self.A[neigh]

      rule = self.__get_interaction_by_types(tv,tn)
      radii = self.__get_particles_radii(tv,tn)

      if nrm < vrad:
        # always reject
        s = 1.0 - nrm/vrad
      elif nrm < radii:
        # attract if rule == 1, else reject
        s = rule*-1.0*\
          (nrm - vrad)/(radii-vrad)
      elif nrm <= vfar:
        # attract if rule == 1, else reject
        s = rule*-1.0*\
          (1.0 - (nrm-radii)/(vfar-radii))
      else:
        continue

      dx = dst[k4]
      dy = dst[k4+1]
      dz = dst[k4+2]

      resx += dx*s
      resy += dy*s
      resz += dz*s

    diffx[v] += resx*stp
    diffy[v] += resy*stp
    diffz[v] += resz*stp

    return num

  @cython.wraparound(False)
  @cython.boundscheck(False)
  @cython.nonecheck(False)
  cpdef long init_rules(
    self,
    long num_particle_types,
    np.ndarray[double, mode="c",ndim=2] interactions,
    np.ndarray[double, mode="c",ndim=2] distances
  ):

    cdef long i
    cdef long j
    cdef long nn = num_particle_types*num_particle_types

    self.interactions = <double *>malloc(nn*sizeof(double))
    self.distances = <double *>malloc(2*num_particle_types*sizeof(double))
    self.num_particle_types = num_particle_types

    for i in xrange(num_particle_types):
      for j in xrange(num_particle_types):
        self.interactions[i*num_particle_types+j] = interactions[i,j]

    for i in xrange(num_particle_types):
      self.distances[i*2] = distances[i,0]
      self.distances[i*2+1] = distances[i,1]

    return 1


  @cython.wraparound(False)
  @cython.boundscheck(False)
  @cython.nonecheck(False)
  cpdef long spawn(
    self,
    np.ndarray[long, mode="c",ndim=1] rnd
  ):

    from numpy.random import random

    cdef long i
    cdef long v
    cdef long a
    cdef long n = len(rnd)

    #with nogil:
    if True:

      for i in xrange(n):

        v = rnd[i]
        self.X[self.vnum] = self.X[v]
        self.Y[self.vnum] = self.Y[v]
        self.Z[self.vnum] = self.Z[v]
        a = self.A[v]
        if a == 0:
          self.A[self.vnum] = 1
        else:
          self.A[self.vnum] = 0

        self.zonemap.__add_vertex(self.vnum)
        self.vnum += 1

    return n


  @cython.wraparound(False)
  @cython.boundscheck(False)
  @cython.nonecheck(False)
  @cython.cdivision(True)
  cpdef long optimize_position(
    self,
    double reject_stp,
    double attract_stp,
    double stp_limit
  ):

    cdef long v
    cdef long i

    cdef double x
    cdef double y
    cdef double z
    cdef double dx
    cdef double dy
    cdef double dz
    cdef double nrm

    cdef long *vertices
    cdef long num
    cdef double *dst

    cdef long asize = self.zonemap.__get_max_sphere_count()

    with nogil, parallel(num_threads=self.procs):
    #with nogil:

      vertices = <long *>malloc(asize*sizeof(long))
      dst = <double *>malloc(asize*sizeof(double)*4)

      for v in prange(self.vnum, schedule='guided'):
      #for v in xrange(self.vnum):

        self.DX[v] = 0.0
        self.DY[v] = 0.0
        self.DZ[v] = 0.0

        num = self.zonemap.__sphere_vertices_dst(
          self.X[v],
          self.Y[v],
          self.Z[v],
          self.__get_particle_far(v),
          vertices,
          dst
        )
        self.__rules(
          v,
          self.DX,
          self.DY,
          self.DZ,
          reject_stp,
          vertices,
          dst,
          num
        )

      free(vertices)
      free(dst)

      for v in prange(self.vnum, schedule='static'):
      #for v in xrange(self.vnum):

        dx = self.DX[v]
        dy = self.DY[v]
        dz = self.DZ[v]

        nrm = sqrt(dx*dx+dy*dy+dz*dz)

        if nrm>stp_limit:
          dx = dx / nrm * stp_limit
          dy = dy / nrm * stp_limit
          dz = dz / nrm * stp_limit

        x = self.X[v] + dx
        y = self.Y[v] + dy
        z = self.Z[v] + dz

        self.X[v] = x
        self.Y[v] = y
        self.Z[v] = z

    with nogil:
      for v in xrange(self.vnum):
        self.zonemap.__update_v(v)

    return 1

  @cython.wraparound(False)
  @cython.boundscheck(False)
  @cython.nonecheck(False)
  cpdef long position_noise(self, np.ndarray[double, mode="c",ndim=2] a):

    cdef long v

    for v in xrange(self.vnum):

      self.X[v] += a[v,0]
      self.Y[v] += a[v,1]
      self.Z[v] += a[v,2]

    return 1

