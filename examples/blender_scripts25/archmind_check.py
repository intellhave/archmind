#  Archmind Non-manifold Geometric Kernel
#  Copyright (C) 2010 Athanasiadis Theodoros
#
#  This software is provided 'as-is', without any express or implied
#  warranty.  In no event will the authors be held liable for any damages
#  arising from the use of this software.
#
#  Permission is granted to anyone to use this software for any purpose,
#  including commercial applications, and to alter it and redistribute it
#  freely, subject to the following restrictions:
#
#  1. The origin of this software must not be misrepresented; you must not
#     claim that you wrote the original software. If you use this software
#     in a product, an acknowledgment in the product documentation would be
#     appreciated but is not required.
#  2. Altered source versions must be plainly marked as such, and must not be
#     misrepresented as being the original software.
#  3. This notice may not be removed or altered from any source distribution.

bl_info = {
    "name": "Mesh Check",
    "author": "Athanasiadis Theodoros",
    "version": (1, 0),
    "blender": (2, 5, 5),
    "api": 31236,
    "location": "View3D > Specials > Mesh Check",
    "description": "Check the mesh",
    "warning": "",
    "wiki_url": "",
    "tracker_url": "",
    "category": "Mesh"}

from archmind_utils import *

def check(askew, amax, dfilename, pfilename):
    mymesh = blender_to_mesh()

    #make a list for min angle triangle distribution
    tri_angles = [0 for i in range(0,12)]

    for f in list(mymesh.faces):
        fa = angles(f)

        try:
            index = int(min( fa ) / (pi * (4.99999/180.0)) )
        except ValueError:
            index = 0

        tri_angles[index] = tri_angles[index] + 1

    if dfilename != '':
        out = open(dfilename, 'w')

        out.write('#Min angle distribution, every 5 degrees\n')
        for i in range(0,12):
            out.write('%3d %d\n' % ((i+1)*5, tri_angles[i]) )

        out.close()

    if pfilename != '':
        out = open(pfilename, 'w')
    
        #for f in mymesh.faces:
        #    for v in f.verts:
        #        out.write('%f %f %f\n' % (v.point.x, v.point.y,v.point.z) )
        #    out.write('\n')

        #for e in mymesh.edges:
        #    out.write('%f %f %f\n' % (e.v0.point.x, e.v0.point.y,e.v0.point.z) )
        #    out.write('%f %f %f\n\n\n' % (e.v1.point.x, e.v1.point.y,e.v1.point.z) )

        for e in mymesh.edges:
            if e.v0.point.y >= -0.05:
                out.write('%f %f %f\n' % (-e.v0.point.x, e.v0.point.y, -e.v0.point.z) )
                #out.write('%f %f %f\n' % (e.v0.point.x, -e.v0.point.y, e.v0.point.z) )
    
            if e.v1.point.y >= -0.05:
                out.write('%f %f %f\n' % (-e.v1.point.x, e.v1.point.y, -e.v1.point.z) )
                #out.write('%f %f %f\n' % (e.v1.point.x, -e.v1.point.y, e.v1.point.z) )

            out.write('\n\n')
    
        out.close()

    return mesh_check(mymesh,askew,askew,amax)

#gui
from bpy.props import *

class CheckClass(bpy.types.Operator):
    bl_idname = 'mesh.check'
    bl_label = 'Check Mesh'
    bl_options = {'REGISTER', 'UNDO'}

    skewangle = bpy.props.FloatProperty(name="Skew Angle :", min=0.0,max=60.0, default=5.0)
    maxangle = bpy.props.FloatProperty(name="Max Angle :", min=60.0,max=180.0, default=170.0)
    anglefile = bpy.props.StringProperty(name="Angle File :", default = "")
    meshfile = bpy.props.StringProperty(name="Mesh File :", default = "")
    
    def execute(self, context):
        msg = str( check(radians( self.skewangle ), radians( self.maxangle ), self.anglefile, self.meshfile) ) 
        self.report( {'INFO'}, msg )
        return {'FINISHED'}

    def invoke(self, context, event):
        wm = context.window_manager
        return wm.invoke_props_dialog(self)

def check_button(self,context):
    self.layout.operator(CheckClass.bl_idname, "Check Mesh")

def register():
    bpy.utils.register_module(__name__)
    bpy.types.VIEW3D_PT_tools_objectmode.prepend(check_button) 

def unregister():
    bpy.utils.unregister_module(__name__)
    bpy.types.VIEW3D_PT_tools_objectmode.remove(check_button) 

if __name__=='__main__':
    register()

