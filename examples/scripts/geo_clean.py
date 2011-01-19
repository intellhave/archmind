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

#Geometry cleaning
#A wavefront of OFF file is expected as an input

import sys
import getopt
from collections import deque
from archmind.geometry import *
from archmind.io import *

def fix_orientation(m):
    """Fixed the faces orientation and returns the number of fixed faces"""
    flips = 0
    checked = {}
    checked = [0 for i in range(0,m.faces_size)]    #build a list with zeros
    queue = deque( [iter(m.faces).next()] )     #add the first face to the queue

    while queue:
        f = queue.popleft()
        checked[f.id] = 1       #mark as fixed

        for ff in f.faces:
            #skip allready fixed faces
            if checked[ff.id]:
                continue

            #get the common edge of the two faces
            e = set(f.edges).intersection(ff.edges).pop()

            if f.edge_cw(e) == ff.edge_cw(e):
                m.flip_face(ff)       #flip face orientation
                flips+=1

            queue.append(ff)

    return flips

def check_duplicate(f):
    """Returns true if the face contains duplicate vertices"""
    return len(set(f.verts)) != f.verts_size


def clean(argv):
    """Removes degenerated faces and fixes faces orientation"""

    try:
        opts, args = getopt.getopt(argv, 's', ['source='])
    except getopt.GetoptError, err:
        print 'exception'
        print str(err)
        sys.exit(1)

    for opt, arg in opts:
        if opt in ('-s', '--source'):
            mesh_filename = arg 
            print 'File = %s' % mesh_filename

    mymesh = mesh()

    print 'Loading mesh... %s' %  mesh_filename
    if not load_from_file(mesh_filename,mymesh):
        sys.exit(1)

    print 'Checking...'

    rem = 0
    for f in list(mymesh.faces):
        if check_duplicate(f):
            mymesh.remove_face(f)
            rem += 1

    #fixing orientation
    flips = fix_orientation(mymesh)

    print 'Faces removed : %d' % rem
    print 'Faces flipped : %d' % flips

    save_to_file(mesh_filename,mymesh)

if __name__ == "__main__":
    clean(sys.argv[1:])
