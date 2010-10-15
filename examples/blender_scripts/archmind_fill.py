#!BPY

"""
Name: 'MeshFill'
Blender: 249
Group: 'Mesh'
Tooltip: 'Fill mesh holes'
"""

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

#Mesh fill tool
#Covers holes of the mesh with perimeter below a certain threshold
#TODO: no check for non planar holes
#TODO: no check for the shape of the hole 
#TODO: at the moment holes may not have more that 128 vertices

from archmind_utils import *

def find(dict,key):
    try:
        return dict[key]
    except KeyError:
        return 0 

def follow(e, checked):
    '''Follows an open boundary'''
    
    vertices = []
    perimeter = 0.0
    cur_edge = e
    start_vert = e.v0
    cur_vert = start_vert
    prev_vert = cur_vert

    while True:
        found_edge = False
        vertices.append( cur_vert )

        #mark the edge as checked
        checked[cur_edge] = True
        
        #calculate the perimeter of the boundary
        perimeter += distance(cur_edge.v0.point,cur_edge.v1.point)
        prev_vert = cur_vert    

        #go to the next point of the edge
        if prev_vert == cur_edge.v0:
            cur_vert = cur_edge.v1
        else:
            cur_vert = cur_edge.v0
        
        #find next edge to follow
        for ee in cur_vert.edges:
            if is_free( ee ) and ee != cur_edge:
            
                #found closed boundary
                if ee == e:
                    return [vertices,perimeter]
            
                #do not select checked edges
                if find(checked,ee):
                    continue

                found_edge = True
                cur_edge = ee
                break

        if not found_edge:
            return [[],0.0]     #failed to find closed boundary


def fill(d):
    '''Fill holes with a perimeter less than d'''
    m = blender_to_mesh()

    checked = {}
    holes = []

    #find the holes of the mesh that satisfy the threshold
    for e in m.edges:
        if is_free(e) and not find(checked,e):
            [hole,perimeter] = follow(e,checked)

            if len(hole) > 0 and perimeter <= d:
                holes.append( hole )

    #cover the holes
    for hole in holes:
        #build a new non convex polygon
        m.add_face( poly(hole) )

    mesh_to_blender(m)
    
#gui
def popup():
    EVT_EXIT = 0

    #Default values
    d = Draw.Create(1000.0)
    
    block = []
    block.append( ('Perimeter', d, 0.0, 1000.0, 'Perimeter' ) )

    if not Draw.PupBlock('Mesh Fill', block):
        return [False,d.val]

    return [True,d.val]

if __name__=='__main__':
    [ok,d] = popup()

    if ok:
        fill(d)

