#!/usr/bin/python
# -*- coding: utf-8 -*-

from __future__ import print_function


def get_args():

  import argparse

  parser = argparse.ArgumentParser()
  parser.add_argument(
    '--procs',
    type=int,
    default=4,
    help='number of processors.'
  )
  parser.add_argument(
    '--nearl',
    type=float,
    default=0.003
  )
  parser.add_argument(
    '--midl',
    type=float,
    default=0.008
  )
  parser.add_argument(
    '--farl',
    type=float,
    default=0.05
  )
  parser.add_argument(
    '--stp',
    type=float,
    default=1.0e-7
  )
  parser.add_argument(
    '--reject',
    type=float,
    default=1.0
  )
  parser.add_argument(
    '--attract',
    type=float,
    default=0.3
  )
  parser.add_argument(
    '--nmax',
    type=int,
    default=1000000
  )
  parser.add_argument(
    '--itt',
    type=int,
    default=10000000000
  )
  parser.add_argument(
    '--vnum',
    type=int,
    default=10000000000
  )
  parser.add_argument(
    '--stat',
    type=int,
    default=100
  )
  parser.add_argument(
    '--export',
    type=int,
    default=1000
  )
  parser.add_argument(
    '--out',
    type=str,
    default='./res/res'
  )
  parser.add_argument(
    '--startRad',
    type=float,
    default=0.01
  )
  parser.add_argument(
    '--startNum',
    type=int,
    default=100
  )

  return parser.parse_args()

def make_info_str(args):
  s = ''
  for k in vars(args):
    s += '# ' + str(k) + ': ' + str(getattr(args,k)) + '\n'
  return s


def print_stats(steps,dm, meta=False):

  from time import strftime
  from time import time

  if isinstance(meta, str):
    meta = ' | {:s}'.format(meta)
  else:
    meta = ''

  print(
    '{:s} | stp: {:d} sec: {:.2f} v: {:d}{:s}'
      .format(
      strftime('%d/%m/%y %H:%M:%S'),
      steps,
      time()-dm.get_start_time(),
      dm.get_vnum(),
      meta
    )
  )

  return

