#!/usr/bin/python
# -*- coding: utf-8 -*-

from __future__ import print_function


def main(args):

  from differentialCloud import DifferentialCloud
  from modules.helpers import print_stats
  from modules.helpers import make_info_str
  from modules.utils import get_initial_cloud
  from modules.utils import export_obj
  from numpy import array

  reject = args.reject*args.stp
  attract = args.attract*args.stp
  stat = args.stat
  export = args.export
  out = args.out
  vnum_max = args.vnum
  start_num = args.startNum
  start_rad = args.startRad

  DC = DifferentialCloud(
    nmax = args.nmax,
    zonewidth = args.farl,
    nearl = args.nearl,
    farl = args.farl,
    procs = args.procs
  )

  xyz, mode = get_initial_cloud(start_num, start_rad)
  rules = array(
    [[1,0],
     [1,0]],
    'int'
  )

  DC.init_cloud(xyz, mode)
  DC.init_rules(rules)

  for i in xrange(args.itt):

    try:

      DC.optimize_position(
        reject,
        attract
      )

      if i%stat==0:
        print_stats(i, DC, meta=None)

      if i%export==0:
        export_obj(
          DC,
          'cloud',
          '{:s}_{:012d}.obj'.format(out, i),
          meta=make_info_str(args)
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

