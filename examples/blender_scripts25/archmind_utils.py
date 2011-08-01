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

from time import clock
from math import *
from archmind.geometry import *
from archmind.io import *
import bpy
import mathutils
import os

def blender_to_mesh():
    '''Extract the blender selected object and returns a mesh'''
    mymesh = mesh()

    sce = bpy.context.scene
    obj = bpy.context.active_object

    blendermesh = obj.to_mesh(bpy.context.scene, True, "PREVIEW")

    vertices = []

    for vert in blendermesh.vertices:
        v = vert.co * obj.matrix_world
        vertices.append( vertex(v[0],v[1],v[2]) )

    for face in blendermesh.faces:
        vert_list = []

        for index in face.vertices:
            vert_list.append( vertices[ index ] )

        mymesh.add_face( poly( vert_list ) )

    return mymesh

# Add the mesh to blender
# edit : Replace existing mesh data
def mesh_to_blender(mymesh,name='ArchMindObject', edit = True):
    '''Convert a mesh to a blender object and adds it to the scene'''
    scene = bpy.context.scene
    obj_act = scene.objects.active

    if edit and not obj_act:
        return
    
    vertices = []
    for v in mymesh.verts:
        vertices.append( [v.point.x,v.point.y,v.point.z ] )
    
    faces = []
    for f in mymesh.faces:
        face_list = []
        for v in f.verts:
            face_list.append( v.id )
   
        #blender does not handle polygons, so we triangulate anything above quad
        if len(face_list) > 4:
            tris = triangulate( f )
            for i in range(0,len(tris),3):
                faces.append( [tris[i].id,tris[i+1].id,tris[i+2].id] )
        elif len(face_list) > 2:
            faces.append( face_list )
    
    newmesh = bpy.data.meshes.new(name)
    newmesh.from_pydata(vertices,[],faces)
    newmesh.update()

    if edit:
        obj = obj_act

        # Replace geometry of existing object
        if obj_act.mode == 'OBJECT':
            old_mesh = obj.data
            obj.data = None
            old_mesh.user_clear()

            # Remove old mesh datablock if no users are left
            if old_mesh.users == 0:
                bpy.data.meshes.remove(old_mesh)

            obj.data = newmesh
            obj.matrix_world = mathutils.Matrix()
    else:
        # Create new object
        obj = bpy.data.objects.new(name, newmesh)
        scene.objects.link(obj)

    obj.select = True

def write_2Dmesh_plot_file(m,filename):
    '''Can write a suitable file that can be converted to postscript in gnuplot, only for 2D (x,y) meshes'''

    postfile = open(filename,'w')
    
    for e in m.edges:
        postfile.write('%f %f' % (e.v0.point.x,e.v1.point.y))

    postfile.close()

def edge_lens(f, decent=1):
    '''Returns a sorted list in decent or ascend order of tuples [ edge_length, edge ] of the face'''
    lens = []

    for e in f.edges:
        lens.append( [distance( e.v0.point, e.v1.point ), e] )

    #return sorted( lens, reverse=decent, cmp=lambda x,y: cmp([x[0],x[1].id],[y[0],y[1].id]) )
    return sorted( lens, reverse=decent )

isqrt3 = 1.0 / sqrt(3.0)

def quality(f):
    r = circum_radius( f.points )
    if isnan(r) or isinf(r):
        return 0.0

    return (isqrt3 * edge_lens(f,decent=0)[0][0]) / circum_radius( f.points )

def angles(f):
    '''Return the list of face angles''' 
    points = list(f.points)

    E0 = (points[1] - points[0]);
    E1 = (points[2] - points[1]);
    E2 = (points[0] - points[2]);

    LE0 = sqrt(E0.x*E0.x+E0.y*E0.y+E0.z*E0.z)
    LE1 = sqrt(E1.x*E1.x+E1.y*E1.y+E1.z*E1.z)
    LE2 = sqrt(E2.x*E2.x+E2.y*E2.y+E2.z*E2.z)

    #normalize
    E0 /= LE0
    E1 /= LE1
    E2 /= LE2

    if LE0 < 0.000001:
        a0 = 0.0
    else:
        try:
            a0 = acos( dot( -E2, E0 ) )
        except ValueError:
            a0 = 0.0
    
    if LE1 < 0.000001:
        a1 = 0.0
    else:
        try:
            a1 = acos( dot( -E0, E1 ) )
        except ValueError:
            a1 = 0.0

    if LE2 < 0.000001:
        a2 = 0.0
    else:
        try:
            a2 = acos( dot( -E1, E2 ) )
        except ValueError:
            a2 = 0.0

    if abs( sum([a0,a1,a2]) - pi ) > 0.01:
        if LE0 > LE1:
            if LE0 > LE2:
                #LE0 is the max
                a2 = pi
            else:
                #LE2 is the max
                a1 = pi
        else:
            if LE1 > LE2:
                #LE1 is the max
                a0 = pi
            else:
                #LE2 is the max
                a1 = pi

        #print 'Degen triangle : ', [a0,a1,a2], ':', [LE0,LE1,LE2]

    return [a0,a1,a2]

def mesh_check(m, angle_needle, angle_skew, angle_max):
    '''Check the mesh and returns stats''' 
    stats = {}
    stats['mesh_quality'] = 0.0
    stats['faces_total'] = m.faces_size
    stats['faces_needles'] = 0
    stats['faces_skew'] = 0
    stats['faces_caps'] = 0
    stats['faces_degen'] = 0
    stats['edges_total'] = m.edges_size
    stats['edges_free'] = 0
    stats['edges_tjoin'] = 0
    stats['verts_total'] = m.verts_size
    stats['verts_free'] = 0

    #face stats
    for f in m.faces:
        fa = angles(f)

        try:
            index = int(min( fa ) / (pi * (4.99999/180.0)) )
        except ValueError:
            stats['faces_degen'] += 1
            index = 0

        if min(fa) < angle_needle:
            stats['faces_needles'] += 1

        if min(fa) < angle_skew:
            stats['faces_skew'] += 1

        if max(fa) > angle_max:
            stats['faces_caps'] += 1

        stats['mesh_quality'] += quality( f ) / m.faces_size

    #edges stats
    for e in m.edges:
        if is_free(e):
            stats['edges_free'] += 1
        elif is_tjoin(e):
            stats['edges_tjoin'] += 1
    
    for v in m.verts:
        if len( list( v.edges ) ) == 0:
            stats['verts_free'] += 1

    return stats


#binary heap
class bheap:
    def __init__(self):
        self.Array = []
        self.Index = {}

    def __len__(self):
        return len(self.Index)

    def push(self,x):
        pos = len(self.Array)
        self.Index[x] = pos
        self.Array.append(x)
        self._percolate_up(pos)

    def clear(self):
        self.Array = []
        self.Index = {}
   
    def pop(self):
        minitem = self.Array[0]
        if len(self.Array) > 1:
            self.Array[0] = self.Array[-1]
            self.Array.pop()

            #update index
            del self.Index[minitem]
            self.Index[self.Array[0]] = 0
            self._percolate_down(0)
        else:
            self.clear()

        return minitem

    def index(self, key):
        return self.Index[key]
   
    def remove_key(self, pos):
        if pos != 0:
            self.change_key(pos,self.Array[0])

        self.pop()

    def change_key(self, pos, newkey):
        oldkey = self.Array[pos]
        del self.Index[oldkey]
        self.Index[newkey] = pos
        self.Array[pos] = newkey

        if oldkey < newkey:
            self._percolate_down(pos)
        elif newkey < oldkey:
            self._percolate_up(pos)

    def remove(self,key):
        try:
            self.remove_key(self.Index[key])
        except KeyError:
            return

    def change(self,key,newkey):
        try:
            self.change_key(self.Index[key],newkey)
        except KeyError:
            return

    def _swap(self,i,j):
        #swap data
        self.Array[i], self.Array[j] = self.Array[j], self.Array[i]
        self.Index[self.Array[i]] = i
        self.Index[self.Array[j]] = j

    def _percolate_up(self,pos):
        while pos >= 0 and self.Array[pos] < self.Array[int(pos / 2)]:
            self._swap(pos,int(pos / 2))
            pos = int(pos / 2)
    
    def _percolate_down(self,pos):
        tmp = self.Array[pos]
        while ((pos * 2) + 1) < len(self.Array):
            child = (pos * 2) + 1
            #find minimun child
            if child < len(self.Array)-1 and self.Array[child+1] < self.Array[child]:
                child += 1
            
            if self.Array[child] < tmp:
                child_item = self.Array[child]
                self.Array[pos] = child_item
                self.Index[child_item] = pos
            else:
                break

            pos = child

        self.Array[pos] = tmp
        self.Index[tmp] = pos
