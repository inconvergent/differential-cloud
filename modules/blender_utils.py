# -*- coding: utf-8 -*-

import bpy

class Obj(object):

  def __init__(self, fn, obj_name):

    self.fn = fn
    self.obj_name = obj_name
    self.obj = self.__import(fn)

    return

  def __import(self, fn):

    bpy.ops.object.select_all(action='DESELECT')

    bpy.ops.import_scene.obj(
      filepath=fn,
      use_smooth_groups=False,
      use_edges=True,
    )

    obj = bpy.context.selected_objects[0]

    return obj

  def get_vertex_color(self):

    from mathutils import Color

    colors = []

    try:

      with open(self.fn+'.x', 'r', encoding='utf8') as f:

        for l in f:
          if l.startswith('#'):
            continue

          values = l.split()
          if not values:
            continue

          if values[0] == 'c':
            c = [float(v) for v in values[1:]]
            colors.append(c)

    except FileNotFoundError:
      return

    mesh = self.obj.data

    if not mesh.vertex_colors:
      mesh.vertex_colors.new()

    col = mesh.vertex_colors.active

    num = len(colors)

    numv = len(self.obj.data.polygons)

    i = 0

    for poly in self.obj.data.polygons:
      loop = poly.loop_indices
      verts = poly.vertices
      for idx,v in zip(loop,verts):
        col.data[idx].color = Color(colors[v])
        i += 1

    print(num, numv, len(col.data), i)

  def move_rescale(self, pos, scale):

    bpy.ops.object.origin_set(type='GEOMETRY_ORIGIN')

    obj = self.obj

    sx,sy,sz = obj.scale

    sx *= scale
    sy *= scale
    sz *= scale

    obj.scale = ((sx,sy,sz))

  def smooth(self, view_levels=1, render_levels=2):

    bpy.context.scene.objects.active = self.obj

    bpy.ops.object.modifier_add(type='SUBSURF')
    self.obj.modifiers['Subsurf'].levels = view_levels
    self.obj.modifiers['Subsurf'].render_levels = render_levels

    bpy.ops.object.shade_smooth()

  def __set_vis(self, frame, vis=True):

    bpy.context.scene.objects.active = self.obj

    bpy.data.scenes['Scene'].frame_current = frame
    bpy.context.active_object.hide = not vis
    bpy.context.active_object.hide_render = not vis

    bpy.context.active_object.keyframe_insert(
      data_path="hide",
      index=-1,
      frame=frame
    )
    bpy.context.active_object.keyframe_insert(
      data_path="hide_render",
      index=-1,
      frame=frame
    )

  def animate_vis(self, ain, aout):

    self.__set_vis(0, False)
    self.__set_vis(ain, True)
    self.__set_vis(aout, False)

  def apply_mat(self):

    mat = bpy.data.materials["Material"]
    self.obj.data.materials.append(mat)

