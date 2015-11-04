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
    double nearl,
    double midl,
    double farl,
    long procs
  ):

    cloud.Cloud.__init__(self, nmax, zonewidth, procs)

    self.nearl = nearl
    self.midl = midl
    self.farl = farl
    self.num_rules = 0

    print('nearl: {:f}'.format(nearl))
    print('midl: {:f}'.format(midl))
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
  cdef long __get_rule_by_types(self, long i, long j) nogil:

    cdef long k = i*self.num_rules+j
    return self.rules[k]

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
    cdef long rule

    for k in range(num):

      neigh = vertices[k]

      if neigh == v:
        continue

      k4 = k*4
      nrm = dst[k4+3]

      if nrm < 1.e-10 or nrm>self.farl:
        continue

      rule = self.__get_rule_by_types(self.A[v], self.A[neigh])

      if nrm < self.nearl:
        # always reject
        #s = self.farl/nrm-1.0
        s = 1.0 - nrm/self.nearl
      elif nrm < self.midl:
        # attract if rule == 1, else reject
        s = pow(-1.0, <double>rule+1.0)*\
          (nrm - self.nearl)/(self.midl-self.nearl)
      elif nrm <= self.farl:
        # attract if rule == 1, else reject
        s = pow(-1.0, <double>rule+1.0)*\
          (1.0 - (nrm-self.midl)/(self.farl-self.midl))
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
    np.ndarray[long, mode="c",ndim=2] r
  ):

    cdef long i
    cdef long j
    cdef long k
    cdef long s = len(r)
    self.rules = <long *>malloc(s*s*sizeof(long))
    self.num_rules = s

    for i in xrange(s):
      for j in xrange(s):
        k = i*s+j
        self.rules[k] = r[i,j]

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
        #self.X[self.vnum] = 0.5 + (1.0-2*random())*0.01
        #self.Y[self.vnum] = 0.5 + (1.0-2*random())*0.01
        #self.Z[self.vnum] = 0.5 + (1.0-2*random())*0.01
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
    cdef double stp_limit = self.nearl*0.3

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
          self.farl,
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

