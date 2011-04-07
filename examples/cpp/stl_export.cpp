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

//Custom stl exporter example

#include "Geometry/Geometry.h"
#include "Io/Io.h"

using namespace arch::geometry;

typedef mesh<> mesh_t;

//Stl exporter class
template<typename mesh_t>
struct Stl
{
    bool can_read(const std::string &filename)const
    {
        return filename.find(".stl") != std::string::npos ||
               filename.find(".STL") != std::string::npos; 
    }
    bool write(const std::string &filename, mesh_t &m)
    {
        std::ofstream Stream(filename.c_str());

        if( !Stream )
        {
            std::cerr << "Failed to open : " << filename << " for parsing\n";
            return false;
        }

        Stream << "solid Archmind\n";
        for( mesh_t::face_iterator_t f = m.faces_begin(); 
            f != m.faces_end(); ++f )
        {
            Stream << "  facet normal " << (*f)->normal() << "\n";
            Stream << "    outer loop\n";

            for( mesh_t::face_t::vertex_iterator_t v = (*f)->verts_begin();
                v != (*f)->verts_end(); ++v )
                Stream << "      vertex " << (*v)->point() << "\n";

            Stream << "    endloop\n";
            Stream << "  endfacet\n";
        }
        Stream << "endsolid Archmind\n";

        return true;
    }
};

//Register custom stl exporter, classes are expected to implement can_read,write
namespace arch{ namespace io{
template<> struct writer_traits<mesh_t>
{
typedef boost::mpl::vector< Stl<mesh_t> > type; 
};
}}

int main(int argc, char **argv)
{
    mesh_t mymesh;

    if( argc != 3 )
    {
        std::cout << argv[0] << " : " << "[input] [output.stl] \n";
        return 1;
    }

    arch::io::load_from_file(argv[1],mymesh);
    arch::io::save_to_file(argv[2],mymesh);

    return 0;
}
