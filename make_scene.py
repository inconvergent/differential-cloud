import bpy


def main(argv):

  from time import time
  from dddUtils.blender import Obj

  fn = argv[0]

  print('importing: ' + fn)

  t1 = time()

  O = Obj(fn,'a')
  O.move_rescale([-0.5]*3, 100)
  O.spheres()
  O.del_mesh()

  print('\ntime:',time()-t1,'\n\n')

  bpy.ops.wm.save_as_mainfile(filepath='./test.blend')


if __name__ == '__main__':

  import sys
  argv = sys.argv
  argv = argv[argv.index("--") + 1:]
  main(argv)

