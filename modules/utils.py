#!/usr/bin/python
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

def load(fn):

  from codecs import open

  vertices = []

  with open(fn, 'r', encoding='utf8') as f:

    for l in f:
      if l.startswith('#'):
        continue

      values = l.split()
      if not values:
        continue
      if values[0] == 'v':
        vertices.append([float(v) for v in values[1:]])

  return {
    'vertices': vertices
  }

def export_obj(dc, obj_name, fn, meta=False):

  from numpy import zeros
  from codecs import open
  from time import time

  vnum = dc.get_vnum()
  np_verts = zeros((vnum,3),'float')

  runtime = time()-dc.get_start_time()

  dc.np_get_vertices(np_verts)

  print('storing mesh ...')
  print('num vertices: {:d}'.format(vnum))

  with open(fn, 'wb', encoding='utf8') as f:

    if meta:
      f.write('# meta:\n')
      f.write(meta+'\n')

    f.write('# info:\n')

    f.write('# vnum: {:d}\n# runtime: {:f}\n\n'
      .format(vnum, runtime))

    f.write('o {:s}\n'.format(obj_name))

    for v in np_verts[:vnum,:]:
      f.write('v {:f} {:f} {:f}\n'.format(*v))

    print('done.')

