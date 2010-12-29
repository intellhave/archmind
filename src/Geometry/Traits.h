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

#ifndef GEOMETRY_TRAITS_H
#define GEOMETRY_TRAITS_H

#include <boost/shared_ptr.hpp>

namespace arch
{

namespace geometry
{

template<typename Traits> class mesh;
template<typename Traits> class face;
template<typename Traits> class edge;
template<typename Traits> class vertex;

typedef std::size_t uid_t;
const uid_t NO_ID = ~uid_t(0);

template< 
    typename Real, 
    typename Vertex, 
    typename Edge, 
    typename Face, 
    typename Mesh 
        >
struct traits
{
    typedef Real real_t;

    typedef Mesh mesh_t;
    typedef Face face_t;
    typedef Edge edge_t;
    typedef Vertex vertex_t;

    typedef boost::shared_ptr< mesh_t > mesh_ptr_t;
 
    typedef boost::shared_ptr< face_t > face_ptr_t;
    typedef boost::shared_ptr< edge_t > edge_ptr_t;
    typedef boost::shared_ptr< vertex_t > vertex_ptr_t;

    static uid_t CounterID;
};

template< typename Real, typename Vertex, typename Edge, typename Face, typename Mesh > 
uid_t traits<Real,Vertex,Edge,Face,Mesh>::CounterID(0);

template<typename Real,typename Derived>
struct traits_base : 
    traits<
        Real,
        vertex< Derived >,
        edge< Derived >,
        face< Derived >,
        mesh< Derived > 
    >
{
};

template<typename Real = float>
struct default_traits : 
    traits_base< 
        Real,
        default_traits<Real> 
        > 
{
};

}

}

#endif
