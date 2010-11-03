#!BPY

"""
Name: 'HermiteSubSurf'
Blender: 249
Group: 'Mesh'
Tooltip: 'Split mesh edges with Hermite Splines'
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

#Hermite based subdivision of surfaces

from archmind_utils import *

def hermite(p0,p1,n0,n1,t):
    '''Calculates the value of the Cubic Hermite Spline for parameter t'''
    #Hermite direction
    d = p0 - p1
    #Hermite tangents
    m0 = cross( cross( d, n0 ), n0 )
    m1 = cross( cross( d, n1 ), n1 )

    t_2 = t*t
    t_3 = t*t*t
    h00 = 2*t_3 - 3*t_2 + 1
    h10 = t_3 - 2*t_2 + t
    h01 = -2*t_3 + 3*t_2
    h11 = t_3 - t_2
    return p0 * h00 + m0 * h10 + p1 * h01 +m1 * h11

def vertex_normal(v):
    '''Area weighted vertex normal'''
    n = vec3(0.0)
    for f in v.faces:
        n += f.normal * area( f.points )
    return normalize(n)

def subsurf():
    m = blender_to_mesh()

    vertex.normal = vec3(0.0)       #extend vertex class to hold the normal 
    #calculate normals
    for v in m.verts:
        v.normal = vertex_normal(v)

    #split the edges of the mesh
    for e in list(m.edges):
        #middle of the Hermite Spline that is defined by the edge endpoints and the corresponding normals
        pt = hermite(e.v0.point,e.v1.point,e.v0.normal,e.v1.normal,0.5)
        #split the edge without triangulation
        v = m.split_edge(e, 0.5)
        #the new point is moved to the middle of the Hermite Spline
        m.set_point(v,pt)

    #the mesh now contains non planar polygons that needs to be converted to triangles/quads 
    #therefore the faces of the mesh are splitted uniformly according to the 
    #following patterns for triangles and quads
    #   p0            p0----p1
    #   /\            |  |  |
    #  /--\            -----
    # /_\/_\          |__|__|
    #p2    p1         p3    p2
    for f in list(m.faces):
        verts = list(f.verts)
        if len(verts) == 6: #triangle
            m.add_face( triangle(verts[5],verts[0],verts[1]) )     #upper triangle
            m.add_face( triangle(verts[1],verts[2],verts[3]) )     #right triangle
            m.add_face( triangle(verts[3],verts[4],verts[5]) )     #left triangle 
            m.add_face( triangle(verts[1],verts[3],verts[5]) )     #middle triangle
        else: #quad
            c = vertex( centroid(verts) )     #quad centroid
            m.add_face( quad(verts[7],verts[0],verts[1],c) )  
            m.add_face( quad(verts[1],verts[2],verts[3],c) )
            m.add_face( quad(verts[3],verts[4],verts[5],c) )
            m.add_face( quad(verts[5],verts[6],verts[7],c) )

        m.remove_face(f)    #remove the old face
    
    mesh_to_blender(m)

if __name__=='__main__':
    subsurf()
