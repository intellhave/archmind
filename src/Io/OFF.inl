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

template<typename mesh_t>
bool OFF<mesh_t>::can_read(const std::string &filename)const
{
    using namespace std;

    return filename.find(".off") != string::npos ||
           filename.find(".OFF") != string::npos;
}

template<typename mesh_t>
bool OFF<mesh_t>::write(const std::string &filename, mesh_t &mesh)
{
    using namespace std;
   
    ofstream Stream(filename.c_str());

    if( !Stream )
    {
        cerr << "off : failed to open : " << filename << " for writing" << endl;
        return false;
    }

    Stream << "OFF\n";

    Stream << "#Exported from Archmind\n";

    Stream << mesh.vertices_size() << " " << mesh.faces_size() << " " << mesh.edges_size() << "\n";

    for( typename mesh_t::vertex_iterator_t v = mesh.vertices_begin(); v != mesh.vertices_end(); ++v )
        Stream << (*v)->point() << "\n";

    for( typename mesh_t::face_iterator_t f = mesh.faces_begin(); f != mesh.faces_end(); ++f )
    {     
        Stream << (*f)->vertices_size() << " ";
       
        for( typename mesh_t::face_t::vertex_iterator_t fv = (*f)->vertices_begin(); fv != (*f)->vertices_end(); ++fv )
            Stream << " " << (*fv)->get_id();
     
        Stream << "\n";       
    }

    return true;
}
