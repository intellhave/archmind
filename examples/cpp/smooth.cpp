//Laplacian smoothing example

#include <Archmind/Geometry/Geometry.h>
#include <Archmind/Geometry/Algorithms.h>
#include <Archmind/Io/Loader.h>

#include <iostream>
#include <vector>

using namespace std;
using namespace arch::geometry;
using namespace arch::io;

typedef mesh<> mesh_t;
typedef mesh_t::mesh_ptr_t mesh_ptr_t;
typedef mesh_t::vertex_t vertex_t;
typedef mesh_t::vertex_ptr_t vertex_ptr_t;
typedef mesh_t::edge_t edge_t;
typedef mesh_t::edge_ptr_t edge_ptr_t;
typedef mesh_t::face_t face_t;
typedef mesh_t::face_ptr_t face_ptr_t;
typedef mesh_t::point_t point_t;

bool is_locked( vertex_ptr_t v )
{
	//Check if any edge of the vertex is free (boundary) or tjoin
	for( vertex_t::edge_iterator_t e = v->edges_begin();
		e != v->edges_end(); ++e )
	{
		if( is_free( *e ) || is_tjoin( *e ) )
			return true;
	}

	return false;
}

int main(int argc, char **argv)
{
	mesh_t mymesh;

	if( argc != 4 )
	{
		std::cerr << argv[0] << " : " << "infile outfile iterations" << std::endl;
		return 1;
	}

	if( !load_from_file(argv[1],mymesh) )
	{
		std::cerr << "Failed to load mesh from file " << argv[1] << std::endl;
		return 1;
	}

	int iterations = atoi(argv[3]);

	//Simultaneous Laplacian smoothing
	for( int i = 0; i < iterations; ++i )
	{
		for( 
			mesh_t::vertex_iterator_t v = mymesh.vertices_begin();
			v != mymesh.vertices_end(); ++v )
		{
			//Move the vertex to the centoid of its neighbors
			if( !is_locked(*v) )
				mymesh.set_point(*v, centroid( (*v)->vertices() ) );
		}
	}

	if( !save_to_file( argv[2], mymesh ) )
		std::cerr << "Failed to write output file " << argv[1] << std::endl;

	return 0;
}
