# -*- coding: utf-8 -*-
# cython: profile=True

from __future__ import division

cimport cython
cimport cloud

from zonemap3d cimport Zonemap3d

from cython.parallel import parallel, prange

from libc.math cimport sqrt
from libc.math cimport pow

from helpers cimport double_array_init
from helpers cimport long_array_init

import numpy as np
cimport numpy as np


cdef class DifferentialCloud(cloud.Cloud):

  def __init__(
    self,
    long nmax,
    double zonewidth,
    double nearl,
    double midl,
    double farl,
    long procs
  ):

    cloud.Cloud.__init__(self, nmax, zonewidth, procs)

    self.nearl = nearl
    self.midl = midl
    self.farl = farl
    self.ptnum = 0

    print('nearl: {:f}'.format(nearl))
    print('farl: {:f}'.format(farl))

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

    return self.interactions[i*self.ptnum+j]

  cdef double __get_interaction_by_id(self, long i, long j) nogil:

    return self.interactions[self.A[i]*self.ptnum+self.A[j]]

  @cython.wraparound(False)
  @cython.boundscheck(False)
  @cython.nonecheck(False)
  cdef double __get_particle_rad(self, long i) nogil:

    return self.distances[2*self.A[i]]

  @cython.wraparound(False)
  @cython.boundscheck(False)
  @cython.nonecheck(False)
  cdef double __get_particles_radii(self, long i, long j) nogil:

    return self.distances[2*self.A[i]] + self.distances[2*self.A[j]]

  @cython.wraparound(False)
  @cython.boundscheck(False)
  @cython.nonecheck(False)
  cdef double __get_particles_harmonic_radii(self, long i, long j) nogil:

    cdef double a = self.distances[2*self.A[i]]
    cdef double b = self.distances[2*self.A[j]]

    return  a*b/(a+b)

  @cython.wraparound(False)
  @cython.boundscheck(False)
  @cython.nonecheck(False)
  cdef double __get_particle_far(self, long i) nogil:

    return self.distances[2*self.A[i]+1]

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

    cdef double dx
    cdef double dy
    cdef double dz
    cdef double nrm
    cdef double s

    cdef double resx = 0.0
    cdef double resy = 0.0
    cdef double resz = 0.0
    cdef double interaction
    cdef double vrad = self.__get_particle_rad(v)
    cdef double vfar = self.__get_particle_far(v)

    cdef double nvmid

    for k in range(num):

      neigh = vertices[k]

      if neigh == v:
        continue

      k4 = k*4
      nrm = dst[k4+3]

      if nrm < 1.e-10 or nrm>vfar:
        continue

      interaction = self.__get_interaction_by_id(v, neigh)
      nvmid = self.__get_particles_harmonic_radii(v, neigh)*4.0
      # TODO: fix this
      nvmid = 0.008

      # print(self.nearl, nvmid, self.midl, self.farl)

      if nrm < vrad:
        # always reject
        s = 1.0 - nrm/vrad
      elif nrm < nvmid:
        # attract if interaction > 0, else reject
        s = interaction*(nrm - vrad)/(nvmid-vrad)
      elif nrm <= vfar:
        # attract if interaction > 0, else reject
        s = interaction*(1.0 - (nrm-nvmid)/(vfar-nvmid))
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
    np.ndarray[double, mode="c",ndim=2] interactions,
    np.ndarray[double, mode="c",ndim=2] distances
  ):

    cdef long i
    cdef long j
    cdef long s = len(interactions)
    self.interactions = <double *>malloc(s*s*sizeof(double))
    self.distances = <double *>malloc(2*s*sizeof(double))
    self.ptnum = s

    for i in xrange(s):
      for j in xrange(s):
        self.interactions[i*s+j] = interactions[i,j]

    for i in xrange(s):
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

    with nogil:
    # if True:

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
    double limit_stp
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
    # with nogil:
    # if True:

      vertices = <long *>malloc(asize*sizeof(long))
      dst = <double *>malloc(asize*sizeof(double)*4)

      for v in prange(self.vnum, schedule='guided'):
      # for v in xrange(self.vnum):

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
      # for v in xrange(self.vnum):

        dx = self.DX[v]
        dy = self.DY[v]
        dz = self.DZ[v]

        nrm = sqrt(dx*dx+dy*dy+dz*dz)

        if nrm>limit_stp:
          nrm = limit_stp/nrm
          dx = dx * nrm
          dy = dy * nrm
          dz = dz * nrm

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

