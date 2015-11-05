# -*- coding: utf-8 -*-

def random_unit_vec(num, scale):

  from numpy.random import normal
  from numpy.linalg import norm
  from numpy import reshape

  rnd = normal(size=(num,3))
  d = norm(rnd,axis=1)
  rnd[:] /= reshape(d, (num,1))

  return rnd*scale

def get_initial_cloud(num, rad):

  from numpy.random import randint

  xyz = random_unit_vec(num, rad)
  xyz += 0.5
  mode = randint(2,size=num)

  return xyz, mode

