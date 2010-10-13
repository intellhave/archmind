#!BPY

"""
Name: 'MeshSplit'
Blender: 249
Group: 'Object'
Tooltip: 'Split mesh edges'
"""
from archmind_utils import *
from collections import deque
from heapq import *

def split(epsilon,maxsplits):
    print 'Spliting mesh...'
    
    mymesh = blender_to_mesh()
    
    heap = []
    
    for e in mymesh.edges:
        heappush( heap, [1.0 / distance(e.v0.point,e.v1.point), e] )

    me = heap[len(heap)-1][0]

    print 'Min edge :', (1.0 / me) * epsilon

    splits = 0
    while heap:
        [l,e] = heappop( heap )

        if l >= me or splits > maxsplits:
            break

        v = mymesh.split_edge(e,0.5,True)
        
        splits += 1

        for e in v.edges:
            heappush( heap, [1.0 / distance(e.v0.point,e.v1.point), e] )

    mesh_to_blender(mymesh)

#gui
def popup():
    EVT_EXIT = 0

	#Default values
    epsilon = Draw.Create(1.0)
    maxsplits = Draw.Create(10000)
	
    block = []
    block.append( ('Epsilon', epsilon, 0.0, 100.0, 'Epsilon' ) )
    block.append( ('Max Splits', maxsplits, 0.0, 100000, 'Max Splits' ) )

    if not Draw.PupBlock('Mesh Split', block):
        return [0,epsilon.val,maxsplits.val]

    return [1,epsilon.val,maxsplits.val]

if __name__=='__main__':
    [pop,epsilon,maxsplits] = popup()

    if pop:
        split(epsilon,maxsplits)

