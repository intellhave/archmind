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

//Vertex
template<typename Traits>
vertex<Traits>::vertex() : m_ID( NO_ID )
{
    m_UniqueID = Traits::CounterID++;

#ifndef NDEBUG
    //std::clog << "Vertex (" << m_UniqueID << ")" << std::endl;
#endif
}

template<typename Traits>
vertex<Traits>::vertex(const point_t &co ) : Point(co), m_ID( NO_ID )
{
    m_UniqueID = Traits::CounterID++;

#ifndef NDEBUG
    //std::clog << "Vertex (" << m_UniqueID << ")" << std::endl;
#endif
}

template<typename Traits>
vertex<Traits>::vertex(const real_t &x, const real_t &y, const real_t &z) : Point(x,y,z), m_ID( NO_ID )
{
    m_UniqueID = Traits::CounterID++;

#ifndef NDEBUG
    //std::clog << "Vertex (" << m_UniqueID << ")" << std::endl;
#endif
}

template<typename Traits>
vertex<Traits>::~vertex()
{

#ifndef NDEBUG
    //std::clog << "~Vertex (" << m_UniqueID << ")" << std::endl;
#endif
}

template<typename Traits>
typename vertex<Traits>::edge_iterator_t vertex<Traits>::edges_begin()const
{
    return Edges.begin();
}

template<typename Traits>
typename vertex<Traits>::edge_iterator_t vertex<Traits>::edges_end()const
{
    return Edges.end();
}

template<typename Traits>
typename vertex<Traits>::vertex_iterator_t vertex<Traits>::verts_begin()const
{
    return vertex_iterator_t( Edges.begin(), unique_id() );
}

template<typename Traits>
typename vertex<Traits>::vertex_iterator_t vertex<Traits>::verts_end()const
{
    return vertex_iterator_t( Edges.end(), NO_ID );
}

template<typename Traits>
std::size_t vertex<Traits>::edges_size()const
{
    return Edges.size();
}

template<typename Traits>
std::size_t vertex<Traits>::verts_size()const
{
    return Edges.size();
}

template<typename Traits>
typename vertex<Traits>::vertex_range_t vertex<Traits>::verts()const
{
    return vertex_range_t( verts_begin(), verts_end() );
}

template<typename Traits>
typename vertex<Traits>::face_iterator_t vertex<Traits>::faces_begin()const
{
    return face_iterator_t( Edges.begin(), Edges.end() );
}

template<typename Traits>
typename vertex<Traits>::face_iterator_t vertex<Traits>::faces_end()const
{
    return face_iterator_t( Edges.end(), Edges.end() );
}

template<typename Traits>
uid_t vertex<Traits>::id()const
{
    return m_ID;
}

template<typename Traits>
void vertex<Traits>::set_id(const uid_t &id)
{
    m_ID = id;
}

template<typename Traits>
uid_t vertex<Traits>::unique_id()const
{
    return m_UniqueID;
}

//Vertex operators
template< typename Traits >
bool vertex<Traits>::operator==( const vertex<Traits> &other )
{
    return m_UniqueID == other.m_UniqueID;
}

template< typename Traits >
bool vertex<Traits>::operator!=( const vertex<Traits> &other )
{
    return m_UniqueID != other.m_UniqueID;
}

template< typename Traits >
bool vertex<Traits>::operator<( const vertex<Traits> &other )
{
    return m_ID < other.m_ID;
}

template< typename Traits >
typename vertex<Traits>::point_t vertex<Traits>::point()const
{
    return Point;
}

//Edge
template<typename Traits>
edge<Traits>::edge() : m_ID(NO_ID)
{
}

template<typename Traits>
edge<Traits>::edge(typename edge<Traits>::vertex_ptr_t v0,typename edge<Traits>::vertex_ptr_t v1) : m_ID( NO_ID )
{
   Vertices[0] = v0;
   Vertices[1] = v1;

   if( v0->unique_id() > v1->unique_id() )
       std::swap( Vertices[0], Vertices[1] );
    
   m_UniqueID = Traits::CounterID++;
   
#ifndef NDEBUG
   //std::clog << "Edge (" << m_UniqueID << ")" << std::endl;
#endif
}


template<typename Traits>
edge<Traits>::~edge()
{
#ifndef NDEBUG
    //std::clog << "~Edge (" << m_UniqueID << ")" << std::endl;
#endif
}

template<typename Traits>
uid_t edge<Traits>::id()const
{
    return m_ID;
}

template<typename Traits>
void edge<Traits>::set_id(const uid_t &id)
{
    m_ID = id;
}

template<typename Traits>
uid_t edge<Traits>::unique_id()const
{
    return m_UniqueID;
}

template<typename Traits>
typename edge<Traits>::vertex_iterator_t edge<Traits>::verts_begin()const
{
    return Vertices.begin();
}

template<typename Traits>
typename edge<Traits>::vertex_iterator_t edge<Traits>::verts_end()const
{
    return Vertices.end();
}

template<typename Traits>
std::size_t edge<Traits>::verts_size()const
{
    return Vertices.size();
}

template<typename Traits>
typename edge<Traits>::vertex_range_t edge<Traits>::verts()const
{
    return vertex_range_t( verts_begin(), verts_end() );
}

template<typename Traits>
typename edge<Traits>::vertex_ptr_t edge<Traits>::v0()const
{
    return Vertices[0];
}

template<typename Traits>
typename edge<Traits>::vertex_ptr_t edge<Traits>::v1()const
{
    return Vertices[1];
}

template<typename Traits>
typename edge<Traits>::face_range_t edge<Traits>::faces()const
{
	return face_range_t( faces_begin(), faces_end() );
}

template<typename Traits>
typename edge<Traits>::face_iterator_t edge<Traits>::faces_begin()const
{
    return Faces.begin();
}

template<typename Traits>
typename edge<Traits>::face_iterator_t edge<Traits>::faces_end()const
{
    return Faces.end();
}

template<typename Traits>
std::size_t edge<Traits>::faces_size()const
{
    return Faces.size();
}

//Edge operators
template< typename Traits >
bool edge<Traits>::operator==( const edge<Traits> &other )
{
    return m_UniqueID == other.m_UniqueID;
}

template< typename Traits >
bool edge<Traits>::operator!=( const edge<Traits> &other )
{
    return m_UniqueID != other.m_UniqueID;
}

template< typename Traits >
bool edge<Traits>::operator<( const edge<Traits> &other )
{
    return m_ID < other.m_ID;
}

//Face
template<typename Traits>
face<Traits>::face() : m_ID(NO_ID)
{
#ifndef NDEBUG
    //std::clog << "Face (" << m_UniqueID << ")" << std::endl;
#endif
}

template<typename Traits>
face<Traits>::~face() 
{
#ifndef NDEBUG
    //std::clog << "~Face (" << m_UniqueID << ")" << std::endl;
#endif
}

template<typename Traits>
typename face<Traits>::point_iterator_t face<Traits>::points_begin()const
{
    return point_iterator_t( verts_begin() );
}

template<typename Traits>
typename face<Traits>::point_iterator_t face<Traits>::points_end()const
{
    return point_iterator_t( verts_end() );
}

template<typename Traits>
std::size_t face<Traits>::points_size()const
{
    return Edges.size();
}

template<typename Traits>
typename face<Traits>::point_range_t face<Traits>::points()const
{
    return point_range_t( points_begin(), points_end() );
}

template<typename Traits>
typename face<Traits>::vertex_iterator_t face<Traits>::verts_begin()const
{
    return vertex_iterator_t( Edges.begin(), EdgesOrientation);
}

template<typename Traits>
typename face<Traits>::vertex_iterator_t face<Traits>::verts_end()const
{
    return vertex_iterator_t( Edges.end(), EdgesOrientation );
}

template<typename Traits>
std::size_t face<Traits>::verts_size()const
{
    return Edges.size();
}

template<typename Traits>
typename face<Traits>::vertex_range_t face<Traits>::verts()const
{
    return vertex_range_t( verts_begin(), verts_end() );
}

template<typename Traits>
typename face<Traits>::edge_iterator_t face<Traits>::edges_begin()const
{
    return Edges.begin();
}

template<typename Traits>
typename face<Traits>::edge_iterator_t face<Traits>::edges_end()const
{
    return Edges.end();
}

template<typename Traits>
std::size_t face<Traits>::edges_size()const
{
    return Edges.size();
}

template<typename Traits>
typename face<Traits>::face_iterator_t face<Traits>::faces_begin()const
{
    return face_iterator_t( Edges.begin(), Edges.end(), m_UniqueID );
}

template<typename Traits>
typename face<Traits>::face_iterator_t face<Traits>::faces_end()const
{
    return face_iterator_t( Edges.end(), Edges.end(), m_UniqueID );
}

//Face operators
template< typename Traits >
bool face<Traits>::operator==( const face<Traits> &other )
{
    return m_UniqueID == other.m_UniqueID;
}

template< typename Traits >
bool face<Traits>::operator!=( const face<Traits> &other )
{
    return m_UniqueID != other.m_UniqueID;
}

template< typename Traits >
bool face<Traits>::operator<( const face<Traits> &other )
{
    return m_UniqueID < other.m_UniqueID;
}

template<typename Traits>
bool face<Traits>::edge_cw( edge_ptr_t e )
{
    typename edge_array_t::iterator ei;
    unsigned pos = 0;

    for( ei = Edges.begin(); ei != Edges.end(); ++ei, ++pos )
    {
        //Correct orientation
        if( *ei == e ) return EdgesOrientation[pos];
    }

    return false;
}

template<typename Traits>
typename face<Traits>::point_t face<Traits>::normal()const
{

    //return math::normal( 
        //Edges[0]->verts()[EdgesOrientation[0] ? 0 : 1]->point(),
        //Edges[1]->verts()[EdgesOrientation[1] ? 0 : 1]->point(),
        //Edges[2]->verts()[EdgesOrientation[2] ? 0 : 1]->point());

    point_t v0 = Edges[0]->verts()[EdgesOrientation[0]]->point();

    //calculate N = E1 x E0
    point_t n = 
        math::cross( 
        Edges[2]->verts()[EdgesOrientation[2]]->point() - v0, 
        Edges[1]->verts()[EdgesOrientation[1]]->point() - v0 );

    //for n-gons, n > 3 use the other points too
    for( std::size_t i = 2; (i+1) < Edges.size(); ++i )
        n += math::cross( 
            Edges[i+1]->verts()[EdgesOrientation[i+1]]->point() - v0, 
            Edges[i]->verts()[EdgesOrientation[i]]->point() - v0 );

    return normalize( n );
}

template<typename Traits>
bool face<Traits>::edge_ccw( edge_ptr_t e )
{
    return !edge_cw( e );
}

template<typename Traits>
void face<Traits>::set_id(const uid_t &id)
{
    m_ID = id;
}

template<typename Traits>
uid_t face<Traits>::unique_id() const
{ 
    return m_UniqueID; 
}

template<typename Traits>
uid_t face<Traits>::id() const
{
    return m_ID;
}

template<typename Traits>
bool face<Traits>::is_triangle()const
{
	return Edges.size() == 3;
}

template<typename Traits>
bool face<Traits>::is_quad()const
{
	return Edges.size() == 4;
}

template<typename Traits>
bool face<Traits>::is_ngon()const
{
	return Edges.size() > 4;
}

