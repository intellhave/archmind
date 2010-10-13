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

#ifndef GEOMETRY_GEOMETRY_H
#define GEOMETRY_GEOMETRY_H

#include "Traits.h"
#include "Iterators.h"

#include <iostream>
#include <bitset>
#include <cstdlib>
#include <cassert>
#include <vector>
#include <deque>
#include <list>
#include <boost/array.hpp>
#include <boost/range/iterator_range.hpp>
#include "../Math/Vector.h"

namespace geometry
{

//!Non manifold surface mesh class
template<typename Traits = default_traits<> >
class mesh
{
public:
	typedef typename Traits::real_t real_t; 
	typedef typename math::vector3<real_t> point_t;

	typedef typename Traits::vertex_t vertex_t;
	typedef typename Traits::edge_t edge_t;
	typedef typename Traits::face_t face_t;
	
    typedef typename Traits::vertex_ptr_t vertex_ptr_t;
	typedef typename Traits::edge_ptr_t edge_ptr_t;
	typedef typename Traits::face_ptr_t face_ptr_t;

	typedef std::vector< vertex_ptr_t > vertex_array_t;
	typedef std::vector< edge_ptr_t > edge_array_t;
	typedef std::vector< face_ptr_t > face_array_t;

	typedef typename vertex_array_t::const_iterator vertex_iterator_t;
	typedef typename edge_array_t::const_iterator edge_iterator_t;
	typedef typename face_array_t::const_iterator face_iterator_t;

	typedef boost::iterator_range< vertex_iterator_t > vertex_range_t;
	typedef boost::iterator_range< edge_iterator_t > edge_range_t;
	typedef boost::iterator_range< face_iterator_t > face_range_t;

	mesh();
    ~mesh();

	vertex_range_t vertices()const;

	vertex_iterator_t vertices_begin()const;
	vertex_iterator_t vertices_end()const;

	std::size_t vertices_size()const;
	
	edge_range_t edges()const;

	edge_iterator_t edges_begin()const;
	edge_iterator_t edges_end()const;

	std::size_t edges_size()const;

	face_range_t faces()const;

	face_iterator_t faces_begin()const;
	face_iterator_t faces_end()const;

	std::size_t faces_size()const;

	//Basic operators

	/*!
	*\brief change the coordinates of a mesh vertex
	*\param v the mesh vertex
	*\param point the coordinates
	*/
	void set_point( vertex_ptr_t v, const point_t &point );

	/*!
	*\brief add a free vertex to the mesh
	*\param v a free vertex (with no connectivity)
	*\return the index of the vertex in the mesh
	*/
	uid_t add_vertex( vertex_ptr_t v );

	/*!
	*\brief removes a free vertex from the mesh
	*\param v
	*\return true if the removal is sucessful
	*/
	bool remove_vertex( vertex_ptr_t v );

	/*!
	*\brief add a face in the mesh
	*\param f a new face
	*\return the index of the face in the mesh
	*\note the face vertices are also added if they are new
	*/
	uid_t add_face( face_ptr_t f );

	/*!
	*\brief removes a face from the mesh
	*\brief f the face to remove
	*\return true if the removal was sucessful
	*/
	bool remove_face( face_ptr_t f );

	/*!
	*\brief splits a given edge into two and creates a new vertex a the point
	*\param e the edge to split
	*\param t the parametric value of the cut (0-1.0)
	*\return the id of the point created 
	*/
	vertex_ptr_t split_edge( edge_ptr_t e, real_t t = 0.5, bool triangulate = false );

	/*!
	*\brief collapses an edge on one of its vertices
	*\param e the edge that will collapse
	*\param v the vertex of the edge e that the edge is collapsed to
	*\return true if the vertex is valid and the collapse operation is succesful
	*/
	bool join_edge( edge_ptr_t e, vertex_ptr_t v );

	/*!
	*\brief creates an edge from the vertices of a polygon and divides the face in two faces
	*\param f the face to split
	*\param v0 a vertex of face f
	*\param v1 a vertex of face f
	*\return the id of the edge created
	*\note the vertices must not form an edge or a degenerated face will be created
	*/
	edge_ptr_t split_face( face_ptr_t f, vertex_ptr_t v0, vertex_ptr_t v1 );

	/*!
	*\brief joins two faces that share a common edge and fuses them together
	*\param f0 a mesh face
	*\param f1 a mesh face
	*\return a fused face the original must share only one common edge
	*\note the order of argument is important for the final face since it inherits the data of f0 but not of f1
	*/
	face_ptr_t join_face( face_ptr_t f0, face_ptr_t f1 );
	//uid_t join_face( face_ptr_t f0, face_ptr_t f1 );

	/*!
	*\flip face orientation
	*\param a mesh face
	*\return always true at the moment
	*/
	bool flip_face( face_ptr_t f );

private:
	//uid_t add_edge( edge_ptr_t e );
	bool remove_edge( edge_ptr_t e );

	void delete_vertex( vertex_ptr_t v );
	void delete_face( face_ptr_t f );
	void delete_edge( edge_ptr_t e );

	vertex_array_t Vertices;
    edge_array_t Edges;
	face_array_t Faces;
};


template<typename Traits>
class face
{
public:
    typedef typename Traits::real_t real_t;
	typedef typename math::vector3<real_t> point_t;

    typedef typename Traits::vertex_t vertex_t;
	typedef typename Traits::edge_t edge_t;
	typedef typename Traits::face_t face_t;
	
	typedef typename Traits::vertex_ptr_t vertex_ptr_t;
	typedef typename Traits::edge_ptr_t edge_ptr_t;
	typedef typename Traits::face_ptr_t face_ptr_t;

	typedef std::vector< edge_ptr_t > edge_array_t;
	typedef std::bitset<Traits::max_poly_pts> bitset_t;

	typedef typename edge_array_t::const_iterator edge_iterator_t;
    typedef face_vertex_iterator<vertex_ptr_t const,edge_iterator_t,bitset_t> vertex_iterator_t;
	typedef face_face_iterator<face_ptr_t const, typename Traits::edge_t::face_iterator_t, edge_iterator_t > face_iterator_t;
    typedef face_point_iterator<point_t,vertex_iterator_t> point_iterator_t;

	typedef boost::iterator_range< point_iterator_t > point_range_t;
	typedef boost::iterator_range< vertex_iterator_t > vertex_range_t;
	typedef boost::iterator_range< edge_iterator_t > edge_range_t;
	typedef boost::iterator_range< face_iterator_t > face_range_t;

    point_iterator_t points_begin()const;
    point_iterator_t points_end()const;

    point_range_t points()const;
	
    std::size_t points_size()const;

	vertex_range_t vertices()const;
	
	vertex_iterator_t vertices_begin()const;
	vertex_iterator_t vertices_end()const;

	std::size_t vertices_size()const;

	edge_iterator_t edges_begin()const;
	edge_iterator_t edges_end()const;

	std::size_t edges_size()const;

	face_iterator_t faces_begin()const;
	face_iterator_t faces_end()const;

	//!Returns true if the edge orientation is clock-wise 
	bool edge_ccw( edge_ptr_t e );

	//!Returns true if the edge orientation is counter clock-wise
	bool edge_cw( edge_ptr_t e );

	point_t normal()const;

	uid_t unique_id() const;
	uid_t get_id() const;

	face();
    ~face();

	template< typename InputIterator >
	face(const InputIterator &start,const InputIterator &end) : m_ID(NO_ID)
	{
		assign(start,end);

		m_UniqueID = Traits::CounterID++;

#ifndef NDEBUG
        std::clog << "Face (" << m_UniqueID << std::endl;
#endif
	}

	// operators
    bool operator==(const face<Traits> &other);
	bool operator!=(const face<Traits> &other);
	
private:
	template<typename> friend class mesh;

	void set_id(const uid_t &id);

	template< typename InputIterator >
	void assign(const InputIterator &start,const InputIterator &end) 
	{
		EdgesOrientation.reset();
		
		//Build the edges
		std::vector< vertex_ptr_t > verts( start, end );
		Edges = edge_array_t( verts.size() );
		
		for( std::size_t i = 0; i < verts.size(); ++i )
		{
			Edges[i] = edge_ptr_t(new edge_t(verts[i],verts[(i+1)%verts.size()]));

			if( Edges[i]->vertices()[0] == verts[i] )
				EdgesOrientation[i] = true;
		}
	}
	edge_array_t Edges;
	bitset_t EdgesOrientation;

	uid_t m_UniqueID;
	uid_t m_ID;
};

template<typename Traits>
class edge
{
public:
    typedef typename Traits::real_t real_t;
	typedef typename math::vector3<real_t> point_t;
    
    typedef typename Traits::vertex_t vertex_t;
	typedef typename Traits::edge_t edge_t;
	typedef typename Traits::face_t face_t;
	
	typedef typename Traits::vertex_ptr_t vertex_ptr_t;
	typedef typename Traits::edge_ptr_t edge_ptr_t;
	typedef typename Traits::face_ptr_t face_ptr_t;
	
    typedef boost::array< vertex_ptr_t, 2 > vertex_array_t;
	typedef std::vector< face_ptr_t > face_array_t;

	typedef typename vertex_array_t::const_iterator vertex_iterator_t;
	typedef typename face_array_t::const_iterator face_iterator_t;

	typedef boost::iterator_range< vertex_iterator_t > vertex_range_t;
	//typedef boost::iterator_range< edge_iterator_t > edge_range_t;
	typedef boost::iterator_range< face_iterator_t > face_range_t;

	edge();
    edge(vertex_ptr_t v0, vertex_ptr_t v1);
	~edge();

	vertex_range_t vertices()const;

	vertex_ptr_t vertices_front()const;
	vertex_ptr_t vertices_back()const;

	vertex_iterator_t vertices_begin()const;
	vertex_iterator_t vertices_end()const;

	std::size_t vertices_size()const;

	face_iterator_t faces_begin()const;
	face_iterator_t faces_end()const;

	std::size_t faces_size()const;

    uid_t unique_id() const; 	
	uid_t get_id() const;

	// operators
    bool operator==(const edge<Traits> &other);
	bool operator!=(const edge<Traits> &other);

private: 
	template<typename> friend class mesh;

	void set_id(const uid_t &id);

	vertex_array_t Vertices;
	face_array_t Faces;

	uid_t m_UniqueID;
    uid_t m_ID;
};

template<typename Traits>
class vertex
{
public:
    typedef typename Traits::real_t real_t;
	typedef typename math::vector3<real_t> point_t;
    
    typedef typename Traits::vertex_t vertex_t;
	typedef typename Traits::edge_t edge_t;
	typedef typename Traits::face_t face_t;
	
	typedef typename Traits::vertex_ptr_t vertex_ptr_t;
	typedef typename Traits::edge_ptr_t edge_ptr_t;
	typedef typename Traits::face_ptr_t face_ptr_t;
	
	typedef std::vector< edge_ptr_t > edge_array_t;

	//Iterators
	//V->E
	typedef typename edge_array_t::const_iterator edge_iterator_t;
    
	//V->F
	typedef vertex_face_iterator<face_ptr_t const,edge_iterator_t, typename Traits::edge_t::face_iterator_t > face_iterator_t;
   
	//V->V
	typedef vertex_vertex_iterator<vertex_ptr_t const,edge_iterator_t> vertex_iterator_t;
    
	typedef boost::iterator_range< vertex_iterator_t > vertex_range_t;

	vertex();
	~vertex();
	vertex(const math::vector3<real_t> &co);
	vertex(const real_t &x, const real_t &y, const real_t &z);

	vertex_range_t vertices()const;

	vertex_iterator_t vertices_begin()const;
	vertex_iterator_t vertices_end()const;

	std::size_t vertices_size()const;

	edge_iterator_t edges_begin()const;
	edge_iterator_t edges_end()const;

	std::size_t edges_size()const;

    face_iterator_t faces_begin()const;
    face_iterator_t faces_end()const;
    
	uid_t unique_id() const; 	
	uid_t get_id() const;

	point_t point()const;
	
	// operators
    bool operator==(const vertex<Traits> &other);
	bool operator!=(const vertex<Traits> &other);
	
private:
	template<typename> friend class mesh;

	void set_id(const uid_t &id);

	point_t Point;
	edge_array_t Edges;

	uid_t m_UniqueID;
	uid_t m_ID;
};

#include "Geometry.inl"
#include "Mesh.inl"

};

#endif
