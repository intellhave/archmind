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

#ifndef GEOMETRY_ALGORITHMS_H
#define GEOMETRY_ALGORITHMS_H

#include <iterator>

namespace arch
{

namespace geometry
{

//!Returns true if the edge is a free edge 
template<typename EdgePtr>
bool is_free( const EdgePtr &e )
{
    return e->faces_size() <= 1;
}

//!Returns true if the edge is a t-join
template<typename EdgePtr>
bool is_tjoin( const EdgePtr &e )
{
    return e->faces_size() >= 3;
}

//!Returns true if the entity is valid (belongs to a mesh)
template<typename EntityPtr>
bool is_valid( const EntityPtr &e )
{
    return e->get_id() != NO_ID;
}

/*!
\brief Calculates the area of a simple polygon
\param Poly a vector of points
\return the signed area of the polygon
*/
template<typename PointVectorType>
typename PointVectorType::value_type::real_t
area(const PointVectorType &poly)
{
    typedef typename PointVectorType::value_type vec_t;
    typedef typename vec_t::real_t real_t;

    const vec_t &v0 = poly[0];
    
    vec_t rf = math::cross( poly[2] - v0, poly[1] - v0 );
    vec_t n( math::normalize( rf ) );       //plane normal

    for( std::size_t i = 2; (i+1) < poly.size(); ++i )
        rf += math::cross( poly[i+1] - v0, poly[i] - v0 );

    return real_t(0.5) * dot( n, rf );
}

template<typename PointIterator>
typename std::iterator_traits< PointIterator >::value_type::real_t
area(const PointIterator &begin, const PointIterator &end)
{
    typedef typename std::iterator_traits< PointIterator >::value_type vec_t;
    
    std::vector< vec_t > poly( begin, end );
    
    return area( poly );
}

template<typename PointIterator>
typename std::iterator_traits< PointIterator >::value_type::real_t
area(const boost::iterator_range< PointIterator> &range)
{
    return area(range.begin(),range.end());
}

template<typename vec_t>
typename vec_t::value_type triangle_area(
       const vec_t &v0, const vec_t &v1, const vec_t &v2 )
{
    return 0.5 * math::magnitude( 
            math::cross(v2 - v0, v1 - v0) );
}

/*!
\brief Checks a planar polygon for convexity
\param Poly a vector of points
\return true if the points form a convex polygon
*/
template<typename PointVectorType>
bool is_convex(const PointVectorType &poly)
{
    typedef typename PointVectorType::value_type vec_t;
    typedef typename vec_t::real_t real_t;

    if( poly.size() <= 3 ) return true;

    vec_t n = math::cross( poly[2] - poly[0], poly[1] - poly[0] );
    
    for( std::size_t i = 1; i < poly.size(); ++i )
    {
        vec_t cs = math::cross( 
            poly[(i+2)%poly.size()] - poly[i], 
            poly[(i+1)%poly.size()] - poly[i] );

        if( dot( cs, n ) < real_t(0.0) )
            return false;
    }

    return true;
}

template<typename PointIterator>
bool is_convex(const PointIterator &begin, const PointIterator &end)
{
    typedef typename std::iterator_traits< PointIterator >::value_type vec_t;
    
    std::vector< vec_t > poly( begin, end );
    
    return is_convex( poly );
}

//!Calculates the circumradius of a triangle
template<typename PointVectorType>
typename PointVectorType::value_type::real_t
circum_radius(const PointVectorType &poly)
{
    using std::sqrt;
    typedef typename PointVectorType::value_type::real_t real_t;

    assert( poly.size() > 1 );

    real_t a = math::distance( poly[0], poly[1] );

    //for edges
    if( poly.size() == 2 )
        return real_t(0.5) * a;

    real_t b = math::distance( poly[1], poly[2] );
    real_t c = math::distance( poly[2], poly[0] );
    
    return (a*b*c)/sqrt( (a+b+c)*(-a+b+c)*(a-b+c)*(a+b-c) );
}

template<typename PointIterator>
typename std::iterator_traits< PointIterator >::value_type::real_t
circum_radius(const PointIterator &begin, const PointIterator &end)
{
    typedef typename std::iterator_traits< PointIterator >::value_type vec_t;
    
    std::vector< vec_t > poly( begin, end );
    return circum_radius( poly );
}

//!Calculates the circumcenter of a triangle
template<typename PointVectorType>
typename PointVectorType::value_type
circum_center(const PointVectorType &poly)
{
    typedef typename PointVectorType::value_type vec_t;
    typedef typename vec_t::real_t real_t;
    using std::pow;

    real_t e0 = math::distance( poly[0], poly[1] );
    real_t e1 = math::distance( poly[1], poly[2] );
    real_t e2 = math::distance( poly[2], poly[0] );

    real_t den = math::magnitude( math::cross(poly[0] - poly[1], poly[1] - poly[2] ) );
    den = real_t(2.0) * den * den;

    real_t a = (e1 * e1 * math::dot( poly[0] - poly[1], poly[0] - poly[2] )) / den;
    real_t b = (e2* e2 * math::dot( poly[1] - poly[0], poly[1] - poly[2] )) / den;
    real_t c = (e0 * e0 * math::dot( poly[2] - poly[0], poly[2] - poly[1] )) / den;

    return (a * poly[0]) + (b * poly[1]) + (c * poly[2]);
}

//!Calculates the circumcenter of a triangle
template<typename PointIterator>
typename std::iterator_traits< PointIterator >::value_type
circum_center(const PointIterator &begin, const PointIterator &end)
{
    typedef typename std::iterator_traits< PointIterator >::value_type vec_t;
    
    std::vector< vec_t > poly( begin, end );
    return circum_center( poly );
}

/*!
\brief Calulcate centroid of list of vertex pointers
\param begin Start iterator
\param end End iterator
\return the centroid
*/
template<typename Iterator>
typename std::iterator_traits<Iterator>::value_type::element_type::point_t
centroid(const Iterator &begin, const Iterator &end)
{
    typedef typename std::iterator_traits<Iterator>::value_type vertex_ptr_t;
    typedef typename vertex_ptr_t::element_type vertex_t;
    typedef typename vertex_t::point_t vec_t;
    typedef typename vertex_t::real_t real_t;

    vec_t c(0.0);
    std::size_t count = 0;
    for( Iterator i = begin; i != end; i++ )
    {
        c += (*i)->point();
        ++count;
    }

    return c / real_t( count );
}

/*!
\brief Calculate centroid of list of vertex pointers
\param range A range iterator
\return the centroid
\*/
template<typename VertexIterator>
typename std::iterator_traits< VertexIterator >::value_type::element_type::point_t
centroid( const boost::iterator_range< VertexIterator > &range )
{
    return centroid( range.begin(), range.end() );
}

template<typename VertexType, typename Iterator>
bool is_ear(const Iterator &begin, const Iterator &end, const VertexType &v0, const VertexType &v1, const VertexType &v2)
{
    //check if the triangle formed contains all the triangles
    for( Iterator j = begin; j != end; ++j )
    {
        if( !j->second || j->first == v0 || j->first == v1 || j->first == v2 )
                continue;
        
        if( math::point_in_triangle( j->first->point(), v0->point(), v1->point(), v2->point() ) )
            return false;
    }

    return true;
}

/*!
\brief Triangulates a non-convex polygonal face using the ear-clipping algorithm
\param f the face to triangulate
\param triangles a list of vertices to store the triangles 
*/
template<typename FacePtr, typename VertexListType >
void triangulate( const FacePtr &f, VertexListType &triangles )
{
    typedef typename FacePtr::element_type face_t;
    typedef typename face_t::point_t vec_t;
    typedef typename face_t::vertex_ptr_t vertex_ptr_t;
    typedef typename face_t::real_t real_t;
    typedef std::pair< vertex_ptr_t, bool > vertex;
    typedef std::list< vertex > list_t;
    
    list_t poly;

    for( typename face_t::vertex_iterator_t i = f->vertices_begin(); i != f->vertices_end(); ++i )
        poly.push_back( vertex( (*i), false ) );

    std::vector< vec_t > poly_points( f->points_begin(), f->points_end() );

    //calculate a poly normal
    const vec_t &p0 = poly_points[0];
    vec_t n = math::cross( poly_points[2] - p0, poly_points[1] - p0 );
    for( std::size_t i = 2; (i+1) < poly_points.size(); ++i )
        n += math::cross( poly_points[i+1] - p0, poly_points[i] - p0 );

    std::size_t reflex_vertices = 0;
    for( typename list_t::iterator i = poly.begin(); i != poly.end(); ++i )
    {
        typename list_t::iterator prev = i, next = i;

        const vertex_ptr_t &v0 = (prev == poly.begin() ? poly.back().first : (--prev)->first );
        const vertex_ptr_t &v1 = i->first;
        const vertex_ptr_t &v2 = (++next == poly.end() ? poly.front().first : next->first );

        vec_t cs = math::cross( v0->point() - v1->point(), v2->point() - v1->point() );
        if( dot( cs, n ) < -math::traits<real_t>::zero_tol )        //reflex vertex
        {
            (*i).second = true;
            ++reflex_vertices;
        }
    }

    //std::clog << "Reflex vertices : " << reflex_vertices << "\n";
    bool ear_found = true;

    while( poly.size() > 3 && ear_found )
    {
        ear_found = false;

        for( typename list_t::iterator i = poly.begin(); i != poly.end(); ++i )
        {
            typename list_t::iterator prev = i, next = i;

            const vertex_ptr_t &v0 = (prev == poly.begin() ? poly.back().first : (--prev)->first );
            const vertex_ptr_t &v1 = i->first;
            const vertex_ptr_t &v2 = (++next == poly.end() ? poly.front().first : next->first );

            vec_t cs = math::cross( v0->point() - v1->point(), v2->point() - v1->point() );
            
            //if( dot( cs, n ) >= real_t(-0.000001) )       //convex vertex
            if( dot( cs, n ) >= -math::traits<real_t>::zero_tol )       //convex vertex
            {       
                if( is_ear(poly.begin(),poly.end(), v0, v1, v2) )
                {
                    ear_found = true;

                    triangles.push_back(v0);
                    triangles.push_back(v1);
                    triangles.push_back(v2);

                    poly.erase( i );
                    break;
                }

                (*i).second = false;        //convex vertex
            }
        }
    }

    if( poly.size() > 3 )
        std::cerr << "Failed to finish triangulation : " << poly.size() << "\n";

    //add the last 3 vertices as a triangle
    typename list_t::iterator i = poly.begin();
    const vertex_ptr_t &v0 = i->first; 
    const vertex_ptr_t &v1 = (++i)->first;
    const vertex_ptr_t &v2 = (++i)->first;
    
    triangles.push_back(v0);
    triangles.push_back(v1);
    triangles.push_back(v2);
}

}

}

#endif
