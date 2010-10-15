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
WaveFront<mesh_t>::WaveFront() 
{
}

template<typename mesh_t>
bool WaveFront<mesh_t>::can_read(const std::string &filename)const
{
    using namespace std;

    return filename.find(".obj") != string::npos ||
           filename.find(".OBJ") != string::npos;
}

template<typename mesh_t>
bool WaveFront<mesh_t>::read(const std::string &filename, mesh_t &mesh)
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
        cerr << "Can not open : " << filename << " for parsing" << endl;
        return false;
    }

    string Line;

    vector< vertex_ptr_t > mesh_vertices;

    while( !Stream.eof() ) 
    {
        vector<string> tokens;

        getline(Stream, Line);
        split(tokens,Line,is_any_of(" "),token_compress_on);

        if( tokens[0] == "v" )      //Vertex     
        {
            real_t x,y,z;

            try
            {
                x = lexical_cast<double>(tokens[1]);
                y = lexical_cast<double>(tokens[2]);
                z = lexical_cast<double>(tokens[3]);
            
                mesh_vertices.push_back( 
                    vertex_ptr_t( new vertex_t(x,y,z) ) );
            }
            catch(bad_lexical_cast &)
            {
                cerr << "WaveFront : Can not safely cast vertex" << endl;
                cerr << "Token : " <<  tokens[1] << "," << tokens[2] << "," << tokens[3] << endl;
            }
        }
        else if( tokens[0] == "vt" )
        {
        }
        else if( tokens[0] == "f" )     //Read Facet
        {
            if( tokens.size() < 4 )
            {
                cerr << "WaveFront : two few face indices!" << endl;
                cerr << "Token : " << Line << endl;
                break;
            }

            std::vector< vertex_ptr_t > vertices;

            try
            {
                vector<string> indextokens;

                //split(tokens,Line,is_any_of(" "),token_compress_on);

                for( std::size_t i = 1; i < tokens.size(); i++ )
                {
                    vector<string> indextokens;

                    //handle uv indices
                    split(indextokens,tokens[i],is_any_of("/"),token_compress_on);

                    std::size_t idx = lexical_cast<int>(indextokens[0])-1;

                    if( idx < mesh_vertices.size() )
                        vertices.push_back( mesh_vertices[idx] );
                    
                    //uv
                    if( indextokens.size() > 1 )
                    {
                        //lexical_cast<int>(indextokens[1])-1   
                    }
                }
            }
            catch(bad_lexical_cast &)
            {
                cerr << "WaveFront : Can not safely cast face indices" << endl;
                break;
            }

            mesh.add_face( 
                face_ptr_t( new face_t( vertices.begin(), vertices.end() ) ) );
        }
    }
    
    return true;
}

template<typename mesh_t>
bool WaveFront<mesh_t>::write(const std::string &filename, mesh_t &mesh)
{
    using namespace std;
    using boost::str;
    using boost::format;

    string Line;
    ofstream Stream(filename.c_str());

    if( !Stream )
    {
        cerr << "wavefront : failed to open : " << filename << " for writing" << endl;
        return false;
    }

    Line = "#Exported from Archmind\n";
    Stream.write( Line.c_str(), Line.size() );
    
    for( typename mesh_t::vertex_iterator_t v = mesh.vertices_begin(); v != mesh.vertices_end(); ++v )
    {
        Line = str( format( "v %f %f %f\n" ) % (*v)->point().x % (*v)->point().y % (*v)->point().z );
        Stream.write( Line.c_str(), Line.size() );
    }

    for( typename mesh_t::face_iterator_t f = mesh.faces_begin(); f != mesh.faces_end(); ++f )
    {
        //write vertices
        Line = "f";
        for( typename mesh_t::face_t::vertex_iterator_t fv = (*f)->vertices_begin(); fv != (*f)->vertices_end(); ++fv )
            Line += str( format(" %d") % ((*fv)->get_id()+1) );

        Line += "\n";
        Stream.write( Line.c_str(), Line.size() );
    }

    return true;
}