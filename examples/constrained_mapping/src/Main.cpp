/*
  Efficient computation of constrained parameterizations on parallel platforms 
  Copyright (C) 2013 Theodoros Athanasiadis, Georgios Zioupos and Ioannis Fudos

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

#include "GeometryTraits.h"
#include "Io/Io.h"
#include "CommandLine.h"
#include "OpenCLSolver.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/foreach.hpp>
#include <deque>

#define foreach BOOST_FOREACH

using namespace parameterization;

bool export_plot_script(const std::string &filename, mesh_t &mesh)
{
     using namespace std;

     ofstream Stream(filename.c_str());
     if( !Stream )
     {
         cerr << "plot : failed to open : " << filename << " for writing" << endl;
         return false;
     }

     Stream << "#Exported from Archmind\n";
  
     foreach( mesh_t::edge_ptr_t e, mesh.edges() ) {
         Stream << e->v0()->u << " " << e->v0()->v << "\n";
         Stream << e->v1()->u << " " << e->v1()->v << "\n\n";
     }

     foreach( mesh_t::vertex_ptr_t v, mesh.verts() ) {
        if( v->pinned ) 
            Stream << v->u << " " << v->v << "\n\n";
     }
}

arch::math::vec3d hsv_to_rgb(double h, double s, double v)
{
    double r, g, b;

    int i = (int)(h * 6.0);
    double f = h * 6. - i;
    double p = v * (1. - s);
    double q = v * (1. - f * s);
    double t = v * (1. - (1. - f) * s);

    switch(i % 6){
        case 0: r = v, g = t, b = p; break;
        case 1: r = q, g = v, b = p; break;
        case 2: r = p, g = v, b = t; break;
        case 3: r = p, g = q, b = v; break;
        case 4: r = t, g = p, b = v; break;
        case 5: r = v, g = p, b = q; break;
    }

    return arch::math::vec3d(r,g,b);
}

arch::math::vec3d red_color( double value )
{
    using std::min;

    double h = 0.0;
    double s = 0.1 + (1.0 - value) * 0.8;
    double v = 1.0;

    return hsv_to_rgb( h, s, v );
}

arch::math::vec3d heatmap_color( double value )
{
    using std::min;

    double h = ((2.0 * value) / 3.0);
    double s = 1.0;
    //double v = 1.0
    double v = 0.75;

    return hsv_to_rgb( h, s, v );
}

bool export_post_script(const std::string &filename,  mesh_t &mesh, bool area_deform = false, int screen_size = 600)
{
     using namespace std;
     using std::min;
     using std::max;
     using namespace arch::math;

     ofstream Stream(filename.c_str());
     if( !Stream )
     {
         cerr << "wavefront : failed to open : " << filename << " for writing" << endl;
         return false;
     }

     //Find bound volume of mesh
     vec2d low(FLT_MAX), up(-FLT_MAX), dim;

     foreach( mesh_t::vertex_ptr_t v, mesh.verts() ) {
         low[0] = min( low[0], v->u );
         low[1] = min( low[1], v->v );

         up[0] = max( up[0], v->u );
         up[1] = max( up[1], v->v );
     }

     dim = up - low;
     double scale_factor = max( dim[0], dim[1] );

     Stream << "%%Exported from Archmind\n";
     Stream << "%%BoundingBox: 0 0 " << screen_size << " " << (screen_size * (dim[1] / dim[0])) << "\n";
  
     double total_area = 0.0;
     double total_det = 0.0;
     foreach( mesh_t::face_ptr_t f, mesh.faces() ) {
         std::vector< vec3d > proj;
         foreach( face_t::vertex_ptr_t fv, f->verts() ) {
             proj.push_back( vec3d( fv->u, fv->v, 0.0 ) );
         }

         vec3d e0 = proj[1] - proj[0];
         vec3d e1 = proj[2] - proj[1];
         vec3d e2 = proj[2] - proj[0];
         double a = fabs( e0.x * e2.y - e0.y * e2.x);
         double A = fabs( area( f->points() ) );
         total_area += A;
         total_det += a;
     }

     foreach( mesh_t::face_ptr_t f, mesh.faces() ) {
         Stream << "newpath\n";
         std::vector< vec3d > proj;
         foreach( face_t::vertex_ptr_t fv, f->verts() ) {
             proj.push_back( vec3d( fv->u, fv->v, 0.0 ) );
         }

         vec3d e0 = proj[1] - proj[0];
         vec3d e1 = proj[2] - proj[1];
         vec3d e2 = proj[2] - proj[0];

         double det = fabs( e0.x * e2.y - e0.y * e2.x );
         double A = fabs( area( f->points() ) );
         
         double px = (proj[0][0] - low[0]) / scale_factor;
         double py = (proj[0][1] - low[1]) / scale_factor;
         Stream << (int)((px * screen_size) + 0.5) << " " << (int)((py * screen_size) + 0.5) << " moveto\n";
         
         for( std::size_t i = 1; i < proj.size(); ++i ) {
             double px = (proj[i][0] - low[0]) / scale_factor;
             double py = (proj[i][1] - low[1]) / scale_factor;
         
             Stream << (int)((px * screen_size) + 0.5) << " " << (int)((py * screen_size) + 0.5) << " lineto\n";
         }

         double area_stretch = fabs( (A/total_area) - (det/total_det) ) / max( fabs(A/total_area), fabs(det/total_det) );
         vec3d rgb = heatmap_color( 1.0 - area_stretch );
         //vec3d rgb = red_color( 1.0 - area_stretch );

         Stream << "closepath\n";
         
         if( !area_deform ) 
            Stream << "0.5 setgray\n";
         else {
            Stream << "gsave\n";
            Stream << rgb[0] << " " << rgb[1] << " " << rgb[2] << " setrgbcolor\n";
            Stream << "fill\n";
            Stream << "grestore\n";
           
         }

         Stream << "0.5 setlinewidth\n";
         Stream << "1 setlinejoin\n";
         Stream << "0.0 setgray\n";
         //Stream << "1.0 0.0 1.0 setrgbcolor\n";
         Stream << "stroke\n";
     }

     foreach( mesh_t::vertex_ptr_t v, mesh.verts() ) {
         if( v->pinned ) { 
             double px = (v->u - low[0]) / scale_factor;
             double py = (v->v - low[1]) / scale_factor;

             Stream << "0 setgray\n";
             Stream << "2.0 setlinewidth\n";
             //Stream << (int)((px * screen_size) + 0.5) << " " << (int)((py * screen_size) + 0.5) << " 4 0 360 arc closepath\n";
             Stream << (int)((px * screen_size) + 0.5) << " " << (int)((py * screen_size) + 0.5) << " 5 0 360 arc closepath\n";
             //Stream << (int)((px * screen_size) + 0.5) << " " << (int)((py * screen_size) + 0.5) << " 6 0 360 arc closepath\n";
             //Stream << (int)((px * screen_size) + 0.5) << " " << (int)((py * screen_size) + 0.5) << " 10 0 360 arc closepath\n";
             Stream << "gsave\n";
             //Stream << "0 0.8 0 setrgbcolor fill\n";
             Stream << "1.0 1.0 0 setrgbcolor fill\n";
             Stream << "grestore\n";
             Stream << "stroke\n";
         }
     }
     
     //Stream << "1 1 0 setrgbcolor fill\n";
     

     //foreach( mesh_t::edge_ptr_t e, mesh.edges() ) {
     //    Stream << e->v0()->u << " " << e->v0()->v << "\n";
     //    Stream << e->v1()->u << " " << e->v1()->v << "\n\n";
     //}
     //
     //foreach( mesh_t::vertex_ptr_t v, mesh.verts() ) {
     //   if( v->pinned ) 
     //       Stream << v->u << " " << v->v << "\n\n";
     //}
}

void read_map_file(const std::string &map_filename, mesh_t &mesh)
{
  using namespace boost;
  using namespace std;

  ifstream Stream(map_filename.c_str());

  if( !Stream ) return;

  string Line;

  int id;
  float u, v;

  while( !Stream.eof() ) 
  {
    vector<string> tokens;

    getline(Stream, Line);	
    trim(Line);
    split(tokens,Line,is_any_of(" "),token_compress_on);

    if( tokens.empty() || tokens.size() != 4 ) continue;

    try
    {
      id  = lexical_cast<int>(tokens[1]);
      u = lexical_cast<float>(tokens[2]);
      v = lexical_cast<float>(tokens[3]);
    } catch(bad_lexical_cast &err)
    {
      cerr << "read_map_file : failed to safely cast tokens : " << err.what() << endl;
      cerr << "tokens : " <<  tokens[1] << "," << tokens[2] << "," << tokens[3] << endl;
      continue;
    }

    if( id >= mesh.verts_size() ) {
        std::cerr << "Wrong input file (id > verts_size)\n";
        return;
    }
    mesh.verts()[id]->pinned = true;
    //mesh.verts()[id]->u = u;
    //mesh.verts()[id]->v = v;
    //std::cout << mesh.verts()[id]->point() << "\n";
    std::cout << "Vertex " << id << " is pinned : ";
    std::cout << u << "," << v << "\n";
    //mesh.set_point( mesh.verts()[id], mesh_t::point_t(u, v,0.0) );
  }
}

//Waverfront with uv exporter class
template<typename mesh_t>
struct WaveFrontUV
{
    bool can_read(const std::string &filename)const
    {
        return filename.find(".obj") != std::string::npos ||
               filename.find(".OBJ") != std::string::npos; 
    }

    void optimize_mesh( mesh_t &mesh )
    {
        std::cout << "Optimizing mesh... ";
        std::vector< bool > added( mesh.faces_size(), false );
        std::vector< mesh_t::face_ptr_t > faces;
        std::vector< mesh_t::vertex_ptr_t > verts( mesh.verts_begin(), mesh.verts_end() );

        for( std::size_t i = 0; i < mesh.faces_size(); )
        {
            std::deque< mesh_t::face_ptr_t > face_queue;
            face_queue.push_front( mesh.faces()[i] );
            added[ i ] = true;
            ++i;
            //std::cout << faces.size() << "\n";
          
            while( !face_queue.empty() ) {
                mesh_t::face_ptr_t f = face_queue.front();
                face_queue.pop_front();
                faces.push_back( f );
             
                foreach( mesh_t::face_t::edge_ptr_t fe, f->edges() ) {
                    foreach( mesh_t::edge_t::face_ptr_t fef, fe->faces() ) {
                        if( !added[ fef->id() ] ) {
                            face_queue.push_back( fef );
                            added[ fef->id() ] = true;
                        }
                    }
                }
            }

            while( added[i] && i < mesh.faces_size() ) ++i;
        }

        //std::cout << "Faces before : " << mesh.faces_size() << "\n";
        //std::cout << "Faces after : " << faces.size() << "\n";

        for( std::size_t i = 0; i < faces.size(); ++i ) 
            mesh.remove_face( faces[i] );

        for( std::size_t i = 0; i < verts.size(); ++i ) 
            mesh.remove_vertex( verts[i] );

        for( std::size_t i = 0; i < faces.size(); ++i ) 
             mesh.add_face( faces[i] );
         std::cout << "Finished\n";
    }

    void optimize_mesh_2( mesh_t &mesh )
    {
        std::cout << "Optimizing mesh... ";
        std::vector< mesh_t::face_ptr_t > faces( mesh.faces_begin(), mesh.faces_end() );
        typedef std::pair< int, mesh_t::vertex_ptr_t > pair_t;
        std::vector< pair_t > verts;

        foreach( mesh_t::vertex_ptr_t v, mesh.verts() ) {
            verts.push_back( pair_t(v->verts_size(), v) );
        }

        std::sort( verts.begin(), verts.end() );

        std::cout << "Faces before : " << mesh.faces_size() << "\n";
        std::cout << "Faces after : " << faces.size() << "\n";

        mesh.clear();

        //for( std::size_t i = 0; i < faces.size(); ++i ) 
        //    mesh.remove_face( faces[i] );

        //for( std::size_t i = 0; i < verts.size(); ++i ) 
        //    mesh.remove_vertex( verts[i].second );
        //
        for( std::size_t i = 0; i < verts.size(); ++i ) 
            mesh.add_vertex( verts[i].second );
        
        for( std::size_t i = 0; i < faces.size(); ++i ) 
             mesh.add_face( faces[i] );

         std::cout << "Finished\n";
    }

    bool write(const std::string &filename, const mesh_t &mesh)
    {
        using namespace std;

        ofstream Stream(filename.c_str());

        if( !Stream )
        {
            cerr << "wavefront : failed to open : " << filename << " for writing" << endl;
            return false;
        }

        Stream << "#Exported from Archmind\n";
        foreach( typename mesh_t::vertex_ptr_t v, mesh.verts() )
            Stream << "v " << v->point() << "\n";
        foreach( typename mesh_t::vertex_ptr_t v, mesh.verts() )
            Stream << "vt " << v->u << " " << v->v << "\n";
        for( typename mesh_t::face_iterator_t f = mesh.faces_begin(); f != mesh.faces_end(); ++f )
        {     
            Stream << "f";
            for( typename mesh_t::face_t::vertex_iterator_t fv = (*f)->verts_begin(); fv != (*f)->verts_end(); ++fv )
                Stream << " " << ((*fv)->id()+1) << "/" << ((*fv)->id()+1);
            Stream << "\n";       
        }

        return true;
    }

    bool read(const std::string &filename, mesh_t &mesh)
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
        vector< real_t > mesh_u;
        vector< real_t > mesh_v;

        size_t line_count = 0;

        while( !Stream.eof() ) 
        {
            vector<string> tokens;

            ++line_count;
            getline(Stream, Line);	
            trim(Line);
            split(tokens,Line,is_any_of(" "),token_compress_on);

            if( tokens.empty() ) continue;

            if( tokens[0] == "v" ){
                if( tokens.size() < 4 )
                {
                    cerr << "wavefront : two few floats!" << endl;
                    cerr << "token : " << Line << endl;
                    break;
                }

                real_t x,y,z;
                try {
                    x = real_t(lexical_cast<float>(tokens[1]));
                    y = real_t(lexical_cast<float>(tokens[2]));
                    z = real_t(lexical_cast<float>(tokens[3]));

                    mesh_vertices.push_back( 
                        vertex_ptr_t( new vertex_t(x,y,z) ) );
                }
                catch(bad_lexical_cast &err) {
                    cerr << "wavefront : failed to safely cast vertex : " << err.what() << endl;
                    cerr << "tokens : " <<  tokens[1] << "," << tokens[2] << "," << tokens[3] << endl;
                }
            }
            else if( tokens[0] == "vt" )
            {
                if( tokens.size() < 3 )
                {
                    cerr << "wavefront : two few floats!" << endl;
                    cerr << "token : " << Line << endl;
                    break;
                }

                real_t u,v;
                try {
                    u = real_t(lexical_cast<float>(tokens[1]));
                    v = real_t(lexical_cast<float>(tokens[2]));
                    
                    mesh_u.push_back( u );
                    mesh_v.push_back( v );
                }
                catch(bad_lexical_cast &err) {
                    cerr << "wavefront : failed to safely cast vertex : " << err.what() << endl;
                    cerr << "tokens : " <<  tokens[1] << "," << tokens[2] << endl;
                }
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
                            std::size_t uv_idx = lexical_cast<int>(indextokens[1])-1;
                            mesh_vertices[idx]->u = mesh_u[ uv_idx ];
                            mesh_vertices[idx]->v = mesh_v[ uv_idx ];
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

        //optimize_mesh_2( mesh );

        return true;
    }
};

//Register custom wavefront exporter, classes are expected to implement can_read,write
namespace arch{ namespace io{
template<> struct writer_traits<mesh_t>
{
typedef boost::mpl::vector< WaveFrontUV<mesh_t> > type; 
};
}}

namespace arch{ namespace io{
template<> struct reader_traits<mesh_t>
{
typedef boost::mpl::vector< WaveFrontUV<mesh_t> > type; 
};
}}

int main(int argc, char *argv[])
{
    cmdline_map_t cmd_args;

    if( !parseCommandLine(argc,argv,cmd_args) ) 
        return 1;

    mesh_t input_mesh;
    mesh_t output_mesh;

    if( !arch::io::load_from_file(cmd_args["source"].as<std::string>(), input_mesh) )
    {
        std::cerr << "Failed to open : " << cmd_args["source"].as<std::string>() << "\n";
        return 1;
    }

    
    foreach( mesh_t::vertex_ptr_t v, input_mesh.verts() ) {
       std::swap( v->u, v->v );
    }


    Options options;

    read_map_file(cmd_args["map"].as<std::string>(), input_mesh);

    options.workgroup = cmd_args["workgroup"].as<int>();
    options.device_type = cmd_args["device"].as<int>();
    options.opt_iters = cmd_args["opt_iters"].as<std::size_t>();
    options.un_iters = cmd_args["un_iters"].as<std::size_t>();
    options.target_residual = cmd_args["res"].as<float>();
    options.weights = 1;
    options.projection_type = cmd_args["proj"].as<int>();
    options.free_boundaries = cmd_args["free"].as<int>();
    std::string param_type = cmd_args["type"].as<std::string>();
    options.scale_iters = cmd_args["scale_iters"].as<std::size_t>();

    if( param_type == "isometric" ) {
        options.energy = 0;
        options.theta = 1.0;
    } else if( param_type == "mips" ) {
        options.energy = 0;
        options.theta = 0.0;
    } else {
        options.energy = 1;
        options.theta = 1.0;    
    }

    if( !options.scale_iters ) options.scale_iters = ((unsigned)~0)>>1;

    int export_ps = cmd_args["ps"].as<int>();

    std::cout.precision(16);

    Stats stats;
    try
    {
        if( export_ps == 1 ) export_plot_script(cmd_args["source"].as<std::string>()+".txt", input_mesh); 
        else if( export_ps > 1 ) export_post_script(cmd_args["source"].as<std::string>()+".eps", input_mesh, export_ps > 2, 2000.0);

        SolverCL solver(input_mesh, output_mesh, options);
        solver.solve(stats);

        if( export_ps == 1 ) export_plot_script(cmd_args["target"].as<std::string>()+".txt", output_mesh); 
        else if( export_ps > 1 ) export_post_script(cmd_args["target"].as<std::string>()+".eps", output_mesh, export_ps > 2, 2000.0);
    }
    catch(const cl::Error &err)
    {
        std::cerr << err.what() << "\n";
        return 0;
    }
    catch(const std::runtime_error &err)
    {
        std::cerr << err.what() << "\n";
        return 0;
    }

    std::cout << "Stats\n";
    std::cout << "Time : " << stats.elapsed_ms << " ms\n";
    std::cout << "Iters : " << stats.iterations << "\n";
    std::cout << "Residual : " << stats.residual << "\n";

    arch::io::save_to_file(cmd_args["target"].as<std::string>(), output_mesh);

    return 0;
}
