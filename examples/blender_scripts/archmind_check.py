#!BPY

"""
Name: 'MeshCheck'
Blender: 249
Group: 'Mesh'
Tooltip: 'Mesh check'
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

def check(askew, amax,dfilename, pfilename):
    print 'Checking mesh...'
    mymesh = blender_to_mesh()

    #make a list for min angle triangle distribution
    tri_angles = [0 for i in range(0,12)]

    for f in list(mymesh.faces):
        fa = angles(f)

        try:
            index = int(min( fa ) / (pi * (4.99999/180.0)) )
        except ValueError:
            index = 0

        tri_angles[index] = tri_angles[index] + 1

    out = open(dfilename, 'w')

    out.write('#Min angle distribution, every 5 degrees\n')
    for i in range(0,12):
        out.write('%3d %d\n' % ((i+1)*5, tri_angles[i]) )

    out.close()

    out = open(pfilename, 'w')
    
    #for f in mymesh.faces:
    #    for v in f.verts:
    #        out.write('%f %f %f\n' % (v.point.x, v.point.y,v.point.z) )
    #    out.write('\n')

    #for e in mymesh.edges:
    #    out.write('%f %f %f\n' % (e.v0.point.x, e.v0.point.y,e.v0.point.z) )
    #    out.write('%f %f %f\n\n\n' % (e.v1.point.x, e.v1.point.y,e.v1.point.z) )

    for e in mymesh.edges:
        if e.v0.point.y >= -0.05:
            out.write('%f %f %f\n' % (-e.v0.point.x, e.v0.point.y, -e.v0.point.z) )
            #out.write('%f %f %f\n' % (e.v0.point.x, -e.v0.point.y, e.v0.point.z) )
    
        if e.v1.point.y >= -0.05:
            out.write('%f %f %f\n' % (-e.v1.point.x, e.v1.point.y, -e.v1.point.z) )
            #out.write('%f %f %f\n' % (e.v1.point.x, -e.v1.point.y, e.v1.point.z) )

        out.write('\n\n')

    
    out.close()

    print mesh_check(mymesh,askew,askew,amax)

    #mesh_to_blender(mymesh)

#gui
def popup():
    EVT_EXIT = 0

	#Default values
    distFileCtrl = Draw.Create('quality.txt')
    postFileCtrl = Draw.Create('mesh.txt')
    skewCtrl = Draw.Create(5.0)
    maxCtrl = Draw.Create(170.0)
	
    block = []
    block.append( ('Skew Angle', skewCtrl, 0.0, 60.0, 'Skew angle limit' ) )
    block.append( ('Max Angle', maxCtrl, 60.0, 180.0, 'Max angle limit' ) )
    block.append( ('Angle File : ', distFileCtrl, 0, 160 ) )  
    block.append( ('Mesh File : ', postFileCtrl, 0, 160 ) )  

    if not Draw.PupBlock('Mesh Check', block):
        return [0,skewCtrl.val,maxCtrl.val,distFileCtrl.val,postFileCtrl.val]

    return [1,radians(skewCtrl.val),radians(maxCtrl.val),distFileCtrl.val,postFileCtrl.val]

if __name__=='__main__':
    [pop,askew,amax,afile,pfile] = popup()

    if pop:
        check(askew,amax,afile,pfile)

