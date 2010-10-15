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


