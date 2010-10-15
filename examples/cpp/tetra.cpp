/*
  Archmind Non-manifold Geometric Kernel
  Copyright (C) 2010 Athanasiadis Theodoros

  This software is provided 'as-is', without any express or implied
  warranty.  In no event will the authors be held liable for any damages
  arising from the use of this software.

  Permission is granted to anyone to use this software for any purpose,
  including commercial applications, and to alter it and redistribute it
  freely, subject to the following restrictions:

  1. The origin of this software must not be misrepresented; you must not
     claim that you wrote the original software. If you use this software
     in a product, an acknowledgment in the product documentation would be
     appreciated but is not required.
  2. Altered source versions must be plainly marked as such, and must not be
     misrepresented as being the original software.
  3. This notice may not be removed or altered from any source distribution.
*/

//A simple example showing how to construct a tetraedron and export is as a wavefront file

#include <Archmind/Geometry/Geometry.h>
#include <Archmind/Io/Loader.h>

#include <vector>

using namespace std;
using namespace arch::geometry;
using namespace arch::io;

typedef mesh<> mesh_t;
typedef mesh_t::mesh_ptr_t mesh_ptr_t;
typedef mesh_t::vertex_t vertex_t;
typedef mesh_t::vertex_ptr_t vertex_ptr_t;
typedef mesh_t::face_t face_t;
typedef mesh_t::face_ptr_t face_ptr_t;
typedef mesh_t::point_t point_t;

int main()
{
    mesh_t tetra;
    vector< vertex_ptr_t > poly(3);
    vector< vertex_ptr_t > vertices(4);

    //Calculate vertices
    vertices[0] = vertex_ptr_t( new vertex_t(-1.0, 0.0,  0.0) );
    vertices[1] = vertex_ptr_t( new vertex_t( 1.0, 0.0,  0.0) );
    vertices[2] = vertex_ptr_t( new vertex_t( 0.0, 0.0,  2.0) );
    vertices[3] = vertex_ptr_t( new vertex_t( 0.0, 2.0,  1.0) );

    //Lower triangle
    poly[0] = vertices[0];
    poly[1] = vertices[1];
    poly[2] = vertices[2];

    tetra.add_face(face_ptr_t( new face_t(poly.begin(),poly.end()) ));

    //Triangle 1
    poly[0] = vertices[0];
    poly[1] = vertices[1];
    poly[2] = vertices[3];

    tetra.add_face(face_ptr_t( new face_t(poly.begin(),poly.end()) ));

    //Triangle 2
    poly[0] = vertices[1];
    poly[1] = vertices[2];
    poly[2] = vertices[3];

    tetra.add_face(face_ptr_t( new face_t(poly.begin(),poly.end()) ));

    //Triangle 3
    poly[0] = vertices[2];
    poly[1] = vertices[0];
    poly[2] = vertices[3];

    tetra.add_face(face_ptr_t( new face_t(poly.begin(),poly.end()) ));

    //Export as a wavefront file
    if( save_to_file("output.obj", tetra) )
        std::cout << "Wavefront file exported (output.obj)" << std::endl;
    else
        std::cerr << "Failed to export mesh (output.obj)" << std::endl;

    return 0;
}
