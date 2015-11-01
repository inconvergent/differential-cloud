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

def export_obj(dm, obj_name, fn, write_intensity=False, meta=False):

  from numpy import zeros
  from codecs import open
  from time import time

  vnum = dm.get_vnum()
  fnum = dm.get_fnum()
  henum = dm.get_henum()
  np_verts = zeros((vnum,3),'float')
  np_tris = zeros((fnum,3),'int')

  runtime = time()-dm.get_start_time()

  dm.np_get_vertices(np_verts)
  dm.np_get_triangles_vertices(np_tris)

  intensity = None

  if write_intensity:
    intensity = zeros(vnum,'double')
    dm.get_vertices_intensity(intensity)

  print('storing mesh ...')
  print('num vertices: {:d}, num triangles: {:d}'.format(vnum, fnum))

  with open(fn, 'wb', encoding='utf8') as f:

    if meta:
      f.write('# meta:\n')
      f.write(meta+'\n')

    f.write('# info:\n')

    f.write('# vnum: {:d}\n# henum: {:d}\n# fnum: {:d}\n# runtime: {:f}\n\n'
      .format(vnum, fnum, henum, runtime))

    f.write('o {:s}\n'.format(obj_name))

    for v in np_verts[:vnum,:]:
      f.write('v {:f} {:f} {:f}\n'.format(*v))

    f.write('s off\n')

    for t in np_tris[:fnum,:]:
      t += 1
      f.write('f {:d} {:d} {:d}\n'.format(*t))

  if write_intensity:

    with open(fn+'.x', 'wb', encoding='utf8') as f:

      f.write('o {:s}\n'.format(obj_name))

      for i in intensity[:vnum]:
        f.write('c {:f} {:f} {:f}\n'.format(*[i]*3))

    print('done.')

