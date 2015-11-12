#!/usr/bin/python
# -*- coding: utf-8 -*-

from __future__ import print_function


def main(args):

  from differentialCloud import DifferentialCloud
  from modules.helpers import print_stats
  from modules.helpers import make_info_str
  from modules.utils import get_initial_cloud
  from dddUtils.ioOBJ import export as export_obj
  from numpy.random import random
  from numpy import array, zeros

  reject = args.reject*args.stp
  attract = args.attract*args.stp
  stat = args.stat
  export = args.export
  out = args.out
  vnum_max = args.vnum
  start_num = args.startNum
  start_rad = args.startRad
  nearl = args.nearl
  midl = args.midl
  farl = args.farl

  np_verts = zeros((args.nmax, 3), 'float')

  DC = DifferentialCloud(
    nmax = args.nmax,
    zonewidth = args.farl,
    nearl = nearl,
    midl = midl,
    farl = farl,
    procs = args.procs
  )

  xyz, mode = get_initial_cloud(start_num, start_rad)

  DC.init_cloud(xyz, mode)
  DC.init_rules(
    array(
      [
        [1.0,1.0],
        [1.0,-1.0]
      ],
      'double'
    ),
    array(
      [
        [nearl,farl],
        [nearl,farl]
      ],
      'double'
    )
  )

  for i in xrange(args.itt):

    try:

      rnd = (random(size=DC.get_vnum())<0.001).nonzero()[0]
      if len(rnd)>0:
        DC.spawn(rnd)

      DC.optimize_position(
        reject,
        attract,
        nearl*0.1
      )

      if i%stat==0:
        print_stats(i, DC, meta=None)

      if i%export==0:

        vnum = DC.np_get_vertices(np_verts)
        export_obj(
          obj_name = 'cloud',
          fn = '{:s}_{:012d}.obj'.format(out, i),
          verts = np_verts[:vnum,:],
          tris = None,
          meta = make_info_str(args)
        )

      if DC.get_vnum()>vnum_max:
        return

    except KeyboardInterrupt:

      break


if __name__ == '__main__':

  from modules.helpers import get_args

  args = get_args()
  print(args)
  main(args)

