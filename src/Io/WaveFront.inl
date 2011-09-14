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
        cerr << "wavefront : failed to open : " << filename << " for parsing" << endl;
        return false;
    }

    string Line;

    vector< vertex_ptr_t > mesh_vertices;
	size_t line_count = 0;

    while( !Stream.eof() ) 
    {
        vector<string> tokens;

		++line_count;
        getline(Stream, Line);	
		trim(Line);
        split(tokens,Line,is_any_of(" "),token_compress_on);

		if( tokens.empty() )
		{
			#ifndef NDEBUG
			cerr << "wavefront : skipping line (" << line_count << ") : " << Line << endl;
			#endif

			continue;
		}

        if( tokens[0] == "v" )      //Vertex     
        {
			if( tokens.size() < 4 )
            {
                cerr << "wavefront : two few floats!" << endl;
                cerr << "token : " << Line << endl;
                break;
            }

            real_t x,y,z;

            try
            {
                x = real_t(lexical_cast<float>(tokens[1]));
                y = real_t(lexical_cast<float>(tokens[2]));
                z = real_t(lexical_cast<float>(tokens[3]));
            
                mesh_vertices.push_back( 
                    vertex_ptr_t( new vertex_t(x,y,z) ) );
            }
            catch(bad_lexical_cast &err)
            {
				cerr << "wavefront : failed to safely cast vertex : " << err.what() << endl;
                cerr << "tokens : " <<  tokens[1] << "," << tokens[2] << "," << tokens[3] << endl;
            }
        }
        else if( tokens[0] == "vt" )
        {
        }
        else if( tokens[0] == "f" )     //Read Facet
        {
            if( tokens.size() < 4 )
            {
                cerr << "wavefront : two few face indices!" << endl;
                cerr << "token : " << Line << endl;
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
                    split(indextokens,tokens[i],is_any_of(" /"),token_compress_on);

                    std::size_t idx = lexical_cast<std::size_t>(indextokens[0])-1;

                    if( idx < mesh_vertices.size() )
                        vertices.push_back( mesh_vertices[idx] );
                    
                    //uv
                    if( indextokens.size() > 1 )
                    {
                        //lexical_cast<int>(indextokens[1])-1   
                    }
                }
            }
            catch(bad_lexical_cast &err)
            {
				cerr << "wavefront : failed to safely cast face indices : " << err.what() << endl;
                break;
            }

            mesh.add_face( 
                face_ptr_t( new face_t( vertices.begin(), vertices.end() ) ) );
        }
    }
    
    return true;
}

template<typename mesh_t>
bool WaveFront<mesh_t>::write(const std::string &filename, const mesh_t &mesh)
{
    using namespace std;
   
    ofstream Stream(filename.c_str());

    if( !Stream )
    {
        cerr << "wavefront : failed to open : " << filename << " for writing" << endl;
        return false;
    }

    Stream << "#Exported from Archmind\n";

    for( typename mesh_t::vertex_iterator_t v = mesh.verts_begin(); v != mesh.verts_end(); ++v )
        Stream << "v " << (*v)->point() << "\n";

    for( typename mesh_t::face_iterator_t f = mesh.faces_begin(); f != mesh.faces_end(); ++f )
    {     
        Stream << "f";
       
        for( typename mesh_t::face_t::vertex_iterator_t fv = (*f)->verts_begin(); fv != (*f)->verts_end(); ++fv )
            Stream << " " << ((*fv)->get_id()+1);
     
        Stream << "\n";       
    }

    return true;
}
