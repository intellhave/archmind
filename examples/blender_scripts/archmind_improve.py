#!BPY

"""
Name: 'MeshImprove'
Blender: 249
Group: 'Mesh'
Tooltip: 'Improves mesh quality'
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

from archmind_utils import *
from collections import deque
from heapq import *
import sys

ITERATIONS = 1
SKEW_ANGLE = radians(5.0)
MAX_ANGLE = radians(175.0)
IMPROVE_ANGLE = radians(14.0)
FEATURE_ANGLE = cos( radians(5.0) )   #dihedral angle between faces
FLIP_FEATURE_ANGLE = cos( radians(15.0) )
NODE_EQUIV = 0.001                    #1mm assuming meters as unit
MIN_VOLUME = 0.001                    #1mm^2 

def matrix_zeros(row,col):
    return [0.0 for i in range(0,row*col)]

def matrix_add(a,b): 
    return [(a[i] + b[i]) for i in range(0,len(a))]

def matrix_solve( m, v ):
    inv = [0.0 for i in range(0,9)]

    inv[0] = m[4]*m[8] - m[5]*m[7];
    inv[1] = m[2]*m[7] - m[1]*m[8];
    inv[2] = m[1]*m[5] - m[2]*m[4];
    inv[3] = m[5]*m[6] - m[3]*m[8];
    inv[4] = m[0]*m[8] - m[2]*m[6];
    inv[5] = m[2]*m[3] - m[0]*m[5];
    inv[6] = m[3]*m[7] - m[4]*m[6];
    inv[7] = m[1]*m[6] - m[0]*m[7];
    inv[8] = m[0]*m[4] - m[1]*m[3];

    det = m[0]*inv[0] + m[1]*inv[3] + m[2]*inv[6]
    
    nv = [0.0,0.0,0.0]

    if isnan(det) or abs(det) < 0.00000001:
        return [False,nv]
       
    idet = 1.0 / det
    for i in range(0,9):
        inv[i] *= idet
        
    for i in range(0,3):
        nv[0] += inv[i] * v[i] 
        nv[1] += inv[3+i] * v[i] 
        nv[2] += inv[6+i] * v[i] 
  
    return [True,nv]

def evaluate(q,v,p):
    return p.x*p.x*q[0] + 2.0*p.x*p.y*q[1] + 2.0*p.x*p.z*q[2] + 2.0*p.x*v[0]\
                        + p.y*p.y*q[4]     + 2.0*p.y*p.z*q[5] + 2.0*p.y*v[1]\
                                           + p.z*p.z*q[8]     + 2.0*p.z*v[2]\
                                                              +         v[3]
          
def init_qv(v):
    #calculate the sum of quadratics for each vertex
    v.q = matrix_zeros( 3, 3 )
    v.v = matrix_zeros( 4, 1 )
   
    for f in v.faces:
        n = f.normal

        r = area( f.points )

        [a,b,c] = [n.x,n.y,n.z]
        d = -dot( n, v.point )
           
        v.q = matrix_add( v.q, [r*a*a, r*a*b, r*a*c, r*a*b, r*b*b, r*b*c, r*a*c, r*b*c, r*c*c] )
        v.v = matrix_add( v.v, [r*a*d, r*b*d, r*c*d, r*d*d] )
   
"""
def check( e ):
    faces = list(e.faces)
    
    if len(faces) != 2:
        return 0

    #arbitrary select a normal
    base_norm = faces[0].normal

    #check the deviation of the face normals around each edge vertices
    for i,v in enumerate(e.verts):
        if is_locked(v):
            return 0

        for f in v.faces:
            if adot( f.normal, base_norm ) < FEATURE_ANGLE:
                return 0
    
    return 1
"""

def is_feature(e):
    e_faces = list(e.faces)
    if len(e_faces) != 2:
        return True
    
    f0 = e_faces[0]
    f1 = e_faces[1]

    if dot(f0.normal,f1.normal) < FEATURE_ANGLE:
        return True

    return False

def volume(points):
    '''tetrahedron volume'''
    m = cross(points[0] - points[2], points[1] - points[2])
    return adot( points[3] - points[2], m ) / 6.0

def flip_edge(m,e,depth=10):
    if depth >= sys.getrecursionlimit()-1:
        print 'flip_edge : reached recursion limit'
        return

    if not is_valid(e):
        #print 'flip edge : invalid'
        return 
    
    e_faces = list(e.faces)

    if len(e_faces) != 2:
        return 

    f0 = e_faces[0]
    f1 = e_faces[1]

    n0 = f0.normal
    n1 = f1.normal

    #check coplanar
    #if min( angles(f0) ) > radians(0.1) and min(angles(f1)) > radians(0.1) and dot(n0,n1) < FEATURE_ANGLE:
    #if magnitude(n0) > 0.99 and magnitude(n1) > 0.99 and adot(n0,n1) < FEATURE_ANGLE:
        #print 'flip edge : Not planar'
    #    return
    
    #find different vertices
    verts = list( set(f0.verts) ^ set(f1.verts) )

    e_verts = list(e.verts)

    #actually this should never happen, only if the faces are degenerated
    if len( verts ) != 2:
        return

    points = [verts[0].point, e_verts[0].point, verts[1].point, e_verts[1].point]
    
    #print 'Volume : %f' % volume( points )

    if min( angles(f0) ) > radians(0.1) and min(angles(f1)) > radians(0.1)\
            and (volume( points ) > MIN_VOLUME or dot(n0,n1) < FLIP_FEATURE_ANGLE):
                return
    #if min( angles(f0) ) > radians(0.1) and min(angles(f1)) > radians(0.1) and dot(n0,n1) < FEATURE_ANGLE:
    #    return

    #check if it is convex
    if not is_convex( points ):
        #print 'flip edge : Not convex'
        return 

    #calculate the circumcircle and the circumradius
    c = circum_center( f0.points )
    r = circum_radius( f0.points )

    #check if we violate the delaunay property
    if (distance(verts[0].point,c)+distance(verts[1].point,c)) >= 2*r:
        #print 'flip edge : Delaunay'
        return

    #create a quad from the two triangles
    f = m.join_face(f0,f1)

    #check if join_face failed
    if f is f0:
        print 'WTF!'
        #m.remove_face(f0)
        return

    #the temporary face has all the surounding edges, so we get a copy
    edges = list(f.edges)

    #split along the other diagonal
    m.split_face(f, verts[0], verts[1])

    #propagate the flipping
    for e in edges:
       flip_edge(m,e,depth+1)


def cap_flip_edge(m,e,cap_queue):
    if not is_valid(e):
        return False
    
    e_faces = list(e.faces)

    if len(e_faces) != 2:
        return False

    f0 = e_faces[0]
    f1 = e_faces[1]

    n0 = f0.normal
    n1 = f1.normal

    #check for coplanar
    if min(angles(f0)) > radians(0.1) and min(angles(f1)) > radians(0.1) and dot(n0,n1) < FEATURE_ANGLE:
    #if dot(n0,n0) > 0.999999 and dot(n1,n1) > 0.999999 and adot(n0,n1) < FEATURE_ANGLE:
        #print 'Feature'
        return False
    
    #find different vertices
    verts = list( set(f0.verts) ^ set(f1.verts) )

    e_verts = list(e.verts)

    if len( verts ) != 2:
        return False

    #check if it is convex
    if not is_convex( [verts[0].point, e_verts[0].point, verts[1].point, e_verts[1].point] ):
        #print 'Not convex'
        return False

    #calculate the circumcircle and the circumradius
    c = circum_center( f0.points )
    r = circum_radius( f0.points )

    #check if we violate the delaunay property
    if (distance(verts[0].point,c)+distance(verts[1].point,c)) >= 2*r:
        #print 'Delaunay'
        return False

    #create a quad from the two triangles
    f = m.join_face(f0,f1)

    #check if join_face failed
    if f is f0:
        print 'WTF!'
        m.remove_face(f0)
        return True

    #split along the other diagonal
    ne = m.split_face(f, verts[0], verts[1])

    for f in ne.faces:
        if is_cap(f):
            cap_queue.append(f)

    return True
    

def mesh_split_edge(m, e, t ):
    #gather different edges
    edges = []
    for f in e.faces:
        for ei in f.edges:
            if ei != e:
                edges.append(ei)
  
    le = is_feature(e)

    v = m.split_edge(e,t,True)
    
    if not le:
        m.set_point(v, centroid( v.verts ) )

    for ee in edges:
        flip_edge(m, ee )

    return v

def is_locked(v):
    for e in v.edges:
        if is_free( e ):
            return 1
    return 0

def choose_edge(f):
    el = edge_lens(f)

    if is_feature(el[0][1]):
        if is_feature(el[1][1]):
            return [el[0][1],False]
        else:
            return [el[0][1],True]

    for nf in el[0][1].faces:
        if nf != f:
            #check if it is terminal
            if el[0][1] == edge_lens(nf)[0][1]:
                #if it is terminal check if the second largest is locked
                if is_feature(el[1][1]):
                    return [el[1][1],True]
                else:
                    #return the terminal edge
                    return [el[0][1],True]
            else:
                return choose_edge(nf)

    #print 'Error'
    return [el[0][1],False]

def check_needle(f):
    lens = edge_lens(f)

    #if max( angles(f) ) < MAX_ANGLE and lens[2][0] / lens[1][0] < 0.25:
    if max( angles(f) ) < MAX_ANGLE or lens[2][0] < NODE_EQUIV:
        return [lens[2][1],True]
    else:
        return [lens[2][1],False]

def find_opt(v0,v1):
    """Finds the quadratic error optimized point"""
    q = matrix_add( v0.q, v1.q )
    v = matrix_add( v0.v, v1.v )

    #check if the points are two close
    if distance(v0.point,v1.point) < NODE_EQUIV:
        return (v0.point + v1.point) * 0.5

    [solved,nv] = matrix_solve(q,v)
    if not solved:
        #print 'Singular matrix'
        #try to find the best between the end points and the centroid
       
        cands = []
        cands.append( [evaluate(q,v,v0.point), v0.point] )
        cands.append( [evaluate(q,v,v1.point), v1.point] )

        c = (v0.point + v1.point) * 0.5
        cands.append( [evaluate(q,v,c), c] )
       
        return min( cands )[1] 

        #return v0.point 

    return -vec3(nv[0],nv[1],nv[2])
    
def is_cap(f):
    return max( angles(f) ) > MAX_ANGLE

def improve(m):

    faces = []
    for f in m.faces:
        faces.append(  [quality(f),f] )
    faces.sort()

    for [q,f] in faces:

        if is_valid(f):

            ma = min( angles(f) )

            if ma < IMPROVE_ANGLE:
                
                max_edge = edge_lens(f)[0][1]

                flip_edge(m,max_edge)

                if not is_valid(f):
                    continue

                #remove skiny with LEPP
                fail_safe = 0
                
                while is_valid(f) and fail_safe < 50:
                    [e,found] = choose_edge(f)
                    if found:
                        #print 'split edge : ', e.v0.point, e.v1.point
                        mesh_split_edge(m,e,0.5)
                        #print 'point : ', sv.point
                        #file.write('%s\n' % sv.point)
                        break

                    fail_safe += 1

    caps = deque([])

    for f in m.faces:
        if( is_cap(f) ):
            caps.append( f )

    #print 'Caps ; ', len(caps)
    last_caps = len(caps)
    
    fail_safe = 0

    while( len(caps) > 0 ) and fail_safe < 10000:
        #print 'caps : ' , len(caps)

        fail_safe += 1

        f = caps.popleft()
            
        if is_valid(f):
            ma = max( angles(f) )
            
            if ma > MAX_ANGLE:
                #print 'Splitting largest edge'
                e = edge_lens(f)[0][1]

                #check if we can flip the edge
                if cap_flip_edge(m,e,caps):
                    continue

                for v in f.verts:
                    if v not in e.verts:
                        v2 = v

                #orthogonal projection
                el = distance(e.v0.point,e.v1.point)

                if el < 0.0001:
                    continue

                t = dot( e.v1.point - e.v0.point, v2.point - e.v0.point ) / (el*el)
                
                n = normalize( e.v1.point - e.v0.point )
                p = e.v0.point + (e.v1.point - e.v0.point) * t

                d = dot( p, n )

                #print 'P : ', n, 'D : ', d

                for of in e.faces:

                    if of != f:
                        slice_face( m, of, e, n, d, caps )
                        break


    faces = []

    for f in m.faces:
        ma = min( angles(f) )
        if ma < SKEW_ANGLE:
            heappush( faces, [ma, f] )
  
    for v in m.verts:
        init_qv(v)
    
    print 'Needles : ', len(faces)

    failed = 0
    #remove needles
    while faces:
        [a,f] = heappop( faces )
        
        if is_valid(f):
            ma = min( angles(f) )
            
            if ma < SKEW_ANGLE:

                [e,is_needle] = check_needle(f)

                if is_locked(e.v0) and is_locked(e.v1):
                    continue
                #collapse edge
                elif is_needle:

                    if is_locked(e.v0):
                        opt = e.v0.point
                    elif is_locked(e.v1):
                        opt = e.v1.point
                    else:
                        opt = find_opt(e.v0,e.v1)

                    op0 = e.v0.point
                    op1 = e.v1.point

                    fail = False

                    #calculate the planes
                    ne = []
                    for v in e.verts:
                        for fv in v.faces:
                            for ef in fv.edges:
                                if e.v0 not in ef.verts and e.v1 not in ef.verts:
                                    pn = normalize( cross( fv.normal, (ef.v1.point - ef.v0.point) ) )
                                    ne.append([v.point,pn,dot(pn,ef.v0.point)])

                    
                    for p in ne:
                        #calculate distance from plane
                        s = dot( p[0], p[1] ) - p[2]
                        sn = dot( opt, p[1] ) - p[2]

                        if s*sn <= 0.0:
                            fail = True
                            failed+=1
                            break

                    if not fail:
                        m.set_point(e.v0,opt)
                        m.set_point(e.v1,opt)

                        m.join_edge(e,e.v0 )
                        e.v0.q = matrix_add(e.v0.q, e.v1.q )
                        e.v0.v = matrix_add(e.v0.v, e.v1.v )
   
                        #Add new needles
                        for vf in e.v0.faces:
                            ma = min( angles(vf) )
                            if ma < SKEW_ANGLE:
                                heappush( faces, [ma, vf] )

    print 'Failed collapses : ', failed
    
def slice_face(m, sf, se, n, d, cap_queue,depth=10):
    #print 'slice face (%d)' % sf.id

    if depth >= sys.getrecursionlimit()-1:
        print 'slice_face : reached recursion limit'
        return
    
    a = max( angles(sf) )

    #check if this is a cap
    if a > MAX_ANGLE:  
        c = -1
    else:
        c = 0

    #keep a copy of edges
    edges = sf.edges

    #split the start edge
    t = clip_line_to_plane(se.v0.point, se.v1.point - se.v0.point, n,d)
        
    if t > 0.000001 and t < 0.99999:
        #v = m.split_edge( se, t, True ) 
        v = m.split_edge( se, t, True ) 
    else:
        #print 'No split occured'
        return

    faces = v.faces
    
    #check the angles of the new faces
    for f in faces:
        a = max( angles(f) )

        if a > MAX_ANGLE:
            c += 1
            cap_queue.append( f )
            #print 'produced cap'

    #if no caps were produced we are fine
    if c == 0:
        return

    #print 'Attempting to continue split'
    for e in edges:
        if is_valid(e):
            t = clip_line_to_plane(e.v0.point, e.v1.point - e.v0.point, n,d)

            if adot( normalize(e.v1.point - e.v0.point), n ) < 1.0 + cos( MAX_ANGLE ):
                #avoid splitting when the edge is almost parallel to the plane
                #print 'Split will create new cap!!'
                return

            if t > 0.000001 and t < 0.99999:
                #find face to go
                for f in e.faces:
                    if f not in faces:
                        slice_face(m, f, e, n , d, cap_queue,depth+1)
                        break

def popup():
    EVT_EXIT = 0
    EVT_CENTER = 1
    EVT_TIME = 2
    EVT_FILE = 3

	#Default values
    iterCtrl = Draw.Create( ITERATIONS )
    skewCtrl = Draw.Create( degrees(SKEW_ANGLE) )
    maxCtrl = Draw.Create( degrees(MAX_ANGLE) )
    improveCtrl = Draw.Create( degrees(IMPROVE_ANGLE) )
    featureCtrl = Draw.Create(  degrees(acos( FEATURE_ANGLE ))  )
    nodeCtrl = Draw.Create(  NODE_EQUIV  )
	
    block = []
    block.append( ('Iterations', iterCtrl, 1, 100, 'Algorithm iterations (3-5 are a good choice)' ) )
    block.append( ('Skew Angle', skewCtrl, 0.0, 60.0, 'Skew angle limit' ) )
    block.append( ('Max Angle', maxCtrl, 60.0, 180.0, 'Max angle limit' ) )
    block.append( ('Improve Angle', improveCtrl, 0.0, 30.0, 'Improve angle limit' ) )
    block.append( ('Feature Angle', featureCtrl, 0.0, 90.0, 'Feature angle limit' ) )
    block.append( ('Node Equiv', nodeCtrl, 0.000001, 1.0, 'Node equivalence threshold' ) )

    if not Draw.PupBlock('Mesh Improve', block):
        return [0,ITERATIONS,SKEW_ANGLE,MAX_ANGLE,IMPROVE_ANGLE,FEATURE_ANGLE,NODE_EQUIV]

    return [1,iterCtrl.val,radians(skewCtrl.val),radians(maxCtrl.val),radians(improveCtrl.val),cos(radians(featureCtrl.val)),nodeCtrl.val]

def check_orientation(m):
    for e in m.edges:
        faces = list(e.faces)

        if len(faces) != 2:
            continue

        if faces[0].edge_cw(e) == faces[1].edge_cw(e):
            return False

    return True

def ent_cmp(self,other):
    if self.id < other.id:
        return -1
    elif self.id == other.id:
        return 0
    else:  
        return +1

def main():
    #extend vertices to contain quadratic data
    vertex.q = matrix_zeros( 3, 3 )
    vertex.v = matrix_zeros( 4, 1 )

    face.__cmp__ = ent_cmp
    edge.__cmp__ = ent_cmp
    vertex.__cmp__ = ent_cmp

    mymesh = blender_to_mesh()

    start_time = clock()

    print 'Initial : ', mesh_check(mymesh,SKEW_ANGLE,IMPROVE_ANGLE,MAX_ANGLE)
    for iter in range(0,ITERATIONS):
        #print 'Orientation correct : ', check_orientation(mymesh)
       
        improve(mymesh)
        print 'Iter : ', iter, mesh_check(mymesh,SKEW_ANGLE,IMPROVE_ANGLE,MAX_ANGLE)
      
    
    mesh_to_blender(mymesh)

    end_time = clock()
    print 'Elapsed %f s' % (end_time - start_time)


if __name__=='__main__':
    [ok,ITERATIONS,SKEW_ANGLE,MAX_ANGLE,IMPROVE_ANGLE,FEATURE_ANGLE,NODE_EQUIV] = popup()

    main()

