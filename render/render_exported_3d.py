import bpy


def main(argv):

  from time import time
  from dddUtils.blender import Cloud

  fn = argv[0]
  img_out = argv[1]

  print(img_out)

  try:
    fn_out = argv[2]
  except Exception:
    fn_out = None

  print('importing: ' + fn)

  t1 = time()

  C = Cloud(fn,'a')
  C.move_rescale(set_pivot=[0.5,-0.5,0.5], pos=[0,0,0], scale=100)
  C.spheres(scale=0.0005)

  print('\ntime:',time()-t1,'\n\n')

  bpy.data.scenes["Scene"].render.filepath = img_out
  bpy.ops.render.render(write_still=True)

  if fn_out:
    bpy.ops.wm.save_as_mainfile(filepath=fn_out)


if __name__ == '__main__':

  import sys
  argv = sys.argv

  argv = argv[argv.index("--") + 1:]
  main(argv)

