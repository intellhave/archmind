import Blender
from Blender import Draw
from time import clock
from math import *
from archmind.geometry import *
from archmind.io import *
import bpy
import os

def blender_to_mesh():
    '''Extract the blender selected object and return a mesh'''
    mymesh = mesh()

    sce = bpy.data.scenes.active

    #get the list of selected objects
    obj_list = list(sce.objects.context)

    if len(obj_list) < 1:
        print 'No object selected'
        return mymesh
    
    obj = sce.objects.active

    blendermesh = obj.getData(mesh=1)
    worldMatrix = obj.matrixWorld

    vertices = []

    for vert in blendermesh.verts:
        v = vert.co * worldMatrix
        vertices.append( vertex(v.x,v.y,v.z) )

    for face in blendermesh.faces:
        vert_list = []

        for v in face.v:
            vert_list.append( vertices[ v.index ] )

        mymesh.add_face( poly( vert_list ) )

    return mymesh


def mesh_to_blender(mymesh,hasuv=False):
    '''Convert a mesh to a blender object and adds it to the scene'''
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
    
    scene = bpy.data.scenes.active
    newmesh = bpy.data.meshes.new('ArchMind')
    newmesh.faceUV = hasuv 
    newmesh.verts.extend(vertices)
    newmesh.faces.extend(faces)
    newobj = scene.objects.new(newmesh, 'ArchMindObject')

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

    return sorted( lens, reverse=decent, cmp=lambda x,y: cmp([x[0],x[1].id],[y[0],y[1].id]) )

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
    '''Check the mesh and return stats''' 
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
