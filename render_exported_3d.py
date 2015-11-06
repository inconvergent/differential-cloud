import bpy


def main(argv):

  from time import time
  from dddUtils.blender import Cloud

  fn = argv[0]

  print('importing: ' + fn)

  t1 = time()

  C = Cloud(fn,'a')
  C.move_rescale(set_pivot=[0.5,-0.5,0.5], pos=[0,0,0], scale=100)
  C.spheres(scale=0.01)
  # O.del_mesh()

  print('\ntime:',time()-t1,'\n\n')

  bpy.ops.wm.save_as_mainfile(filepath='./test.blend')


if __name__ == '__main__':

  import sys
  argv = sys.argv
  argv = argv[argv.index("--") + 1:]
  main(argv)

