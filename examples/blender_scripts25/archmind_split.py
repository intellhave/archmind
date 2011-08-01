
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

#Hermite based subdivision of surfaces

bl_info = {
    "name": "UniformSplit",
    "author": "Athanasiadis Theodoros",
    "version": (1, 0),
    "blender": (2, 5, 5),
    "api": 31236,
    "location": "View3D > Specials > Unisplit",
    "description": "Split mesh edges with Hermite Splines",
    "warning": "",
    "wiki_url": "",
    "tracker_url": "",
    "category": "Mesh"}

from archmind_utils import *
from heapq import *

def split_mesh(epsilon,maxsplits):
    mymesh = blender_to_mesh()
    
    heap = []
   
    #build a heap with edge lengths
    for e in mymesh.edges:
        heappush( heap, [1.0 / distance(e.v0.point,e.v1.point), e] )

    #min edge
    me = heap[-1][0]

    print('Min edge :', (1.0 / me) * epsilon)

    splits = 0
    while heap:
        #pop the largest edge
        [l,e] = heappop( heap )

        if l >= me or splits > maxsplits:
            break

        #split the edge and store the vertex created
        v = mymesh.split_edge(e,0.5,True)
        
        splits += 1

        #add the new edges to the heap
        for e in v.edges:
            heappush( heap, [1.0 / distance(e.v0.point,e.v1.point), e] )

    mesh_to_blender(mymesh)

#GUI
from bpy.props import *

class SplitClass(bpy.types.Operator):
    bl_idname = 'mesh.unisplit'
    bl_label = 'Uniform Split'
    bl_options = {'REGISTER', 'UNDO'}

    epsilon = bpy.props.FloatProperty(name='Epsilon', default=1.0, min = 0.0, max = 100.0)
    maxsplits = bpy.props.IntProperty(name='Max Splits', default=1000, min = 0, max = 100000)

    def execute(self, context):
        split_mesh(self.epsilon,self.maxsplits)
        return {'FINISHED'}

def split_button(self, context):
    self.layout.operator(SplitClass.bl_idname, text="Uniform Split")

def register():
    bpy.utils.register_module(__name__)
    bpy.types.VIEW3D_PT_tools_objectmode.append(split_button)

def unregister():
    bpy.utils.unregister_module(__name__)
    bpy.types.VIEW3D_PT_tools_objectmode.remove(split_button)

if __name__=='__main__':
    register()

