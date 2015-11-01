# -*- coding: utf-8 -*-
# cython: profile=True

from __future__ import division

cimport cython
from libc.stdlib cimport malloc, free

from cython.parallel import parallel, prange

from libc.math cimport sqrt
from libc.math cimport pow as cpow
from libc.math cimport cos
from libc.math cimport sin
from libc.math cimport fabs
from libc.math cimport M_PI

from helpers cimport long_array_init
from helpers cimport double_array_init

import numpy as np
cimport numpy as np
cimport cython


cdef double TWOPI = M_PI*2


def dict_list_add(dict d, k, v):

  if d.has_key(k):
    d[k].append(v)
  else:
    d[k] = [v]


cdef class Cloud:
  """
  """

  def __init__(self, long nmax, double zonewidth, long procs):
    """
    """

    from time import time

    self.nmax = nmax
    self.vnum = 0

    self.zonewidth = zonewidth

    self.procs = procs
    self.state = 0

    self.nz = <long>(1.0 /zonewidth)
    self.start_time = time()

    if self.nz<3:
      self.nz = 1
      self.zonewidth = 1.0

    self.zonemap = Zonemap3d(self.nz)
    self.zonemap.__assign_xyz_arrays(self.X, self.Y, self.Z)

    print('nmax: {:d}'.format(nmax))
    print('number of zones: {:d}'.format(self.nz))
    print('zonewidth: {:f}'.format(zonewidth))

    return

  def __cinit__(self, long nmax, *arg, **args):

    self.X = <double *>malloc(nmax*sizeof(double))
    double_array_init(self.X,nmax,0.)

    self.Y = <double *>malloc(nmax*sizeof(double))
    double_array_init(self.Y,nmax,0.)

    self.Z = <double *>malloc(nmax*sizeof(double))
    double_array_init(self.Z,nmax,0.)

    self.A = <long *>malloc(nmax*sizeof(long))
    long_array_init(self.A,nmax,0)

    return

  def __dealloc__(self):

    free(self.X)
    free(self.Y)
    free(self.Z)
    free(self.A)

    return

  @cython.wraparound(False)
  @cython.boundscheck(False)
  @cython.nonecheck(False)
  cdef long __valid_new_vertex(self, double x, double y, double z) nogil:
    """
    check that x,y is within unit square
    """

    if x<0. or x>1.:
      return -1

    if y<0. or y>1.:
      return -1

    if z<0. or z>1.:
      return -1

    return 1

  @cython.wraparound(False)
  @cython.boundscheck(False)
  @cython.nonecheck(False)
  cdef void __set_vertex_type(self, long v1, long a) nogil:

    self.A[v1] = a

    return

  @cython.wraparound(False)
  @cython.boundscheck(False)
  @cython.nonecheck(False)
  cdef long __new_vertex(self, double x, double y, double z) nogil:

    if self.__valid_new_vertex(x,y,z)<0:
      return -1

    cdef long vnum = self.vnum

    self.X[vnum] = x
    self.Y[vnum] = y
    self.Z[vnum] = z

    self.zonemap.__add_vertex(vnum)

    self.vnum += 1
    return vnum

  @cython.wraparound(False)
  @cython.boundscheck(False)
  @cython.nonecheck(False)
  cpdef long np_get_vertices(self, np.ndarray[double, mode="c",ndim=2] x):

    cdef long v
    cdef long c

    for v in xrange(self.vnum):

      if self.A[v]>-1:

        x[c, 0] = self.X[v]
        x[c, 1] = self.Y[v]
        x[c, 2] = self.Z[v]
        c += 1

    return c

  cpdef long set_rules(
    self,
    np.ndarray[double, mode="c",ndim=2] r
  ):

    cdef long i
    self.rules = r

    return 1

  cpdef long init_cloud(
    self,
    np.ndarray[double, mode="c",ndim=2] x,
    np.ndarray[double, mode="c",ndim=1] a
  ):

    cdef long i

    for i in xrange(len(a)):

      self.X[i] = x[i, 0]
      self.Y[i] = x[i, 1]
      self.Z[i] = x[i, 2]
      self.A[i] = a[i]

  @cython.nonecheck(False)
  cpdef long get_vnum(self):

    return self.vnum

