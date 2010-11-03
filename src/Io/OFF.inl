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
bool OFF<mesh_t>::read(const std::string &filename, mesh_t &mesh)
{
    using namespace boost;
    using namespace std;

    typedef typename mesh_t::vertex_ptr_t vertex_ptr_t;
    typedef typename mesh_t::vertex_t vertex_t;

    typedef typename mesh_t::face_ptr_t face_ptr_t;
    typedef typename mesh_t::face_t face_t;

    typedef typename mesh_t::real_t real_t;

    ifstream Stream(filename.c_str());

    if( !Stream )
    {
        cerr << "off : failed to open : " << filename << " for parsing" << endl;
        return false;
    }

    string Line;

    getline(Stream, Line);
    
    if( Line.substr(0,3) != "OFF" )
    {
        cerr << "off : failed to read OFF tag : " << filename << '\n';
        return false;
    }

    std::size_t vertices = 0, faces = 0;

    //skip comments and read vertices and faces
    while( !Stream.eof() ) 
    {
        vector<string> tokens;
        getline(Stream, Line);
		trim(Line);
        split(tokens,Line,is_any_of(" "),token_compress_on);

        //skip comments
        if( !tokens.empty() && tokens[0][0] == '#' )
            continue;

        else if( tokens.size() == 3 )
        {
            try
            {
                vertices = lexical_cast<size_t>(tokens[0]);
                faces = lexical_cast<size_t>(tokens[1]);
            }
            catch(bad_lexical_cast &)
            {
                cerr << "off : failed to safely cast sizes : " << Line << '\n';     
                return false;
            }

            break;
        }
        else
        {
            cerr << "off : failed to read faces and vertices : " << Line << '\n';
            return false;
        }
    }

    while( !Stream.eof() ) 
    {
        vector<string> tokens;
        getline(Stream, Line);
		trim(Line);
        split(tokens,Line,is_any_of(" "),token_compress_on);

        if( vertices > 0 && tokens.size() == 3 )      //Read Vertex     
        {
            real_t x,y,z;

            try
            {
                x = real_t(lexical_cast<float>(tokens[0]));
                y = real_t(lexical_cast<float>(tokens[1]));
                z = real_t(lexical_cast<float>(tokens[2]));
            
                mesh.add_vertex( vertex_ptr_t( new vertex_t(x,y,z) ) );               
            }
            catch(bad_lexical_cast &)
            {
                cerr << "off : failed to safely cast vertex\n";
                cerr << "tokens : " <<  tokens[0] << "," << tokens[1] << "," << tokens[2] << '\n';
            }

            vertices--;
        }
        else if( faces > 0 && !tokens.empty() )    //Read Facet
        {
            std::size_t npts, idx;
            std::vector< vertex_ptr_t > vertices;   //Facet vertices
          
            try
            {
                //Read number of points
                npts = lexical_cast<int>(tokens[0]);

                assert(npts < tokens.size());

                for( std::size_t i = 0; i < npts; i++ )
                {            
                    idx = lexical_cast<int>(tokens[i+1]);               
                    
                    assert( idx < mesh.vertices_size() );
                   
                    vertices.push_back( mesh.vertices()[idx] );                                  
                }
            }
            catch(bad_lexical_cast &)
            {
                cerr << "off : failed to read face indices : " << Line << '\n';
                break;
            }

            mesh.add_face( 
                face_ptr_t( new face_t( vertices.begin(), vertices.end() ) ) );
        
            faces--;
        }
    }
    
    return faces == 0 && vertices == 0;
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
