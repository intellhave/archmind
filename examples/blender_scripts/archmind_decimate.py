#!BPY

"""
Name: 'MeshDecimate'
Blender: 249
Group: 'Object'
Tooltip: 'Quadric-Based Decimation tool'
"""

#Quadric-Based decimation tool
#See 'Quadric-Based Polygonal Surface Simplification' by Michael Garland for the background theory
#TODO: edge contraction does not correctly update the affected edges
#TODO: can_move function is not efficiently implemented since it is reevaluated all the time

from heapq import *
from archmind_utils import *

NODE_EQUIV = 0.001      #node equivalence limit 
DECIMATE_NUM = 1000     #decimation target limit
FEATURE_ANGLE = cos( radians(60.0) )   #dihedral angle between faces that define features

def matrix_zeros(row,col):
    '''Matrix set'''
    return [0.0 for i in range(0,row*col)]

def matrix_add(a,b): 
    '''Matrix addition'''
    return [(a[i] + b[i]) for i in range(0,len(a))]

def matrix_solve( m, v ):
    '''Invert matrix m and multiply with vector v'''
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
    '''Evaluate quadric error'''
    return p.x*p.x*q[0] + 2.0*p.x*p.y*q[1] + 2.0*p.x*p.z*q[2] + 2.0*p.x*v[0]\
                        + p.y*p.y*q[4]     + 2.0*p.y*p.z*q[5] + 2.0*p.y*v[1]\
                                           + p.z*p.z*q[8]     + 2.0*p.z*v[2]\
                                                              +         v[3]
          
def init_qv(v):
    '''Calculate the sum of quadratics for a vertex v'''
    v.q = matrix_zeros( 3, 3 )
    v.v = matrix_zeros( 4, 1 )
   
    for f in v.faces:
        n = f.normal
    
        #weight by the area
        r = area( f.points )

        [a,b,c] = [n.x,n.y,n.z]
        d = -dot( n, v.point )
           
        v.q = matrix_add( v.q, [r*a*a, r*a*b, r*a*c, r*a*b, r*b*b, r*b*c, r*a*c, r*b*c, r*c*c] )
        v.v = matrix_add( v.v, [r*a*d, r*b*d, r*c*d, r*d*d] )
  
def find_opt(v0,v1):
    '''Finds the quadratic error optimized point of the edge v0,v1'''
    q = matrix_add( v0.q, v1.q )
    v = matrix_add( v0.v, v1.v )

    #check if the points are two close
    if distance(v0.point,v1.point) < NODE_EQUIV:
        opt = (v0.point + v1.point) * 0.5
        return [evaluate(q,v,opt),opt]

    [solved,nv] = matrix_solve(q,v)
    if not solved:
        #Singular matrix, try to find the best between the end points and the centroid
       
        cands = []
        cands.append( [evaluate(q,v,v0.point), v0.point] )
        cands.append( [evaluate(q,v,v1.point), v1.point] )
        c = (v0.point + v1.point) * 0.5
        cands.append( [evaluate(q,v,c), c] )
       
        return min( cands )

    opt = -vec3(nv[0],nv[1],nv[2])
    return [evaluate(q,v,opt),opt]
  
def can_move(v):
    '''Returns false if the vertex is a feature'''
    for e in v.edges:
        if is_feature( e ):
            return False
    return True

def is_feature(e):
    '''Returns true if the edge is a feature edge'''
    e_faces = list(e.faces)
    if len(e_faces) != 2:
        return True
    
    f0 = e_faces[0]
    f1 = e_faces[1]

    if dot(f0.normal,f1.normal) < FEATURE_ANGLE:
        return True

    return False

def decimate(m,target):
    '''Decimate a mesh m until a certain face threshold is reached'''

    #calculate quadrics
    for v in m.verts:
        init_qv(v)
   
    edges = []
    #calculate quadrics errors
    for e in m.edges:
        if can_move(e.v0) and can_move(e.v1):
            #find optimized point and add it to the heap
            heappush( edges, [find_opt(e.v0,e.v1)[0], e] )

    while edges and m.faces_size > target:
        #pop the min  
        [err,e] = heappop(edges)
       
        #check if the edge has been invalidated
        if not is_valid(e):
            continue

        #find the optimized point
        opt = find_opt(e.v0,e.v1)[1]
        fail = False

        #calculate the planes around the two points to check for element inversion
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

            #detect element inversion
            if s*sn <= 0.0:
                fail = True
                break

        if not fail:
            #move the points of the edge to the optimized point
            m.set_point(e.v0,opt)
            m.set_point(e.v1,opt)

            #new point
            np = e.v0

            #update quadric matrices
            np.q = matrix_add(e.v0.q, e.v1.q )
            np.v = matrix_add(e.v0.v, e.v1.v )
          
            #collapse the degenerated edge to its first endpoint
            m.join_edge(e,e.v0 )

            #add new edges to the heap
            for ne in np.edges:
                if can_move(ne.v0) and can_move(ne.v1):
                    #find optimized point and add it to the heap
                    heappush( edges, [find_opt(ne.v0,ne.v1)[0], ne] )

def popup():
    EVT_EXIT = 0

	#Default values
    epsilon = Draw.Create(NODE_EQUIV)
    target = Draw.Create(DECIMATE_NUM)
	
    block = []
    block.append( ('Epsilon', epsilon, 0.0, 1.0, 'Epsilon' ) )
    block.append( ('Target faces', target, 0.0, 100000, 'Target faces' ) )

    if not Draw.PupBlock('Mesh Decimate', block):
        return [0,epsilon.val,target.val]

    return [1,epsilon.val,target.val]

def main():
    #extend vertices to contain quadratic data
    vertex.q = matrix_zeros( 3, 3 )
    vertex.v = matrix_zeros( 4, 1 )

    mymesh = blender_to_mesh()
    decimate(mymesh,DECIMATE_NUM)
    mesh_to_blender(mymesh)

if __name__=='__main__':
    [ok,NODE_EQUIV,DECIMATE_NUM] = popup()

    if ok:
        main()

