#Laplacian smoothing example
#A wavefront file is expected as an input

import sys
import getopt
from archmind.geometry import *
from archmind.io import *

try:
	opts, args = getopt.getopt(sys.argv[2:], 's', ['source='])
except getopt.GetoptError, err:
	print str(err)
	sys.exit(1)

for opt, arg in opts:
	if opt in ('-s', '--source'):
		mesh_filename = arg

#helper function to check if a vertex can be moved
def locked(v):
	for e in v.edges:
		if is_free( e ):
			return 1
	return 0

mymesh = mesh()

print 'Loading mesh...'
if not load_from_file(mesh_filename,mymesh):
	sys.exit(1)

print 'Smoothing...'
NumOfIters = 20

for i in range(0,NumOfIters):
	for v in mymesh.verts:
		if not locked(v):
			mymesh.set_point( v, centroid( v.verts ) )

print 'Saving mesh...'
save_to_file(mesh_filename,mymesh)


