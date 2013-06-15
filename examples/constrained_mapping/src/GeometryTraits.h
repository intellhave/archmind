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

#ifndef SPHERE_GEOMETRY_H
#define SPHERE_GEOMETRY_H

#include "Geometry/Geometry.h"
#include "Geometry/Algorithms.h"

namespace parameterization
{
template<typename Traits>
struct param_face : arch::geometry::face<Traits>
{
  typedef typename Traits::real_t real_t;

  param_face() : face() {}
  
  template< typename InputIterator >
  param_face(const InputIterator &start,const InputIterator &end)
  {
    face::face(start,end);
  }

  real_t cosA;
  real_t cosB;
  real_t cosC;
};

template<typename Traits>
struct param_vertex : arch::geometry::vertex<Traits>
{
  typedef typename Traits::point_t point_t;
  typedef typename Traits::real_t real_t;

  param_vertex(const point_t &co) : pinned(false), arch::geometry::vertex<Traits>(co) {}
  param_vertex(const real_t &x, const real_t &y, const real_t &z) : pinned(false), arch::geometry::vertex<Traits>(x,y,z) {}
  point_t normal()const
  {
    point_t n(0.0f);
    for( face_iterator_t f = faces_begin();
         f != faces_end(); ++f )
    {
      n += (*f)->normal();
    }
    return -normalize(n);
  }

  std::vector< float > Weights;
  bool pinned;
  real_t u,v;
};

struct param_traits : 
	arch::geometry::traits<
        float,
        param_vertex< param_traits >,
		    arch::geometry::edge< param_traits >,
        arch::geometry::face< param_traits >,
        arch::geometry::mesh< param_traits > 
    >
{
};

typedef arch::geometry::mesh<param_traits> mesh_t;
typedef mesh_t::point_t point_t;

typedef mesh_t::vertex_t vertex_t;
typedef mesh_t::edge_t edge_t;
typedef mesh_t::face_t face_t;
typedef mesh_t::vertex_ptr_t vertex_ptr_t;
typedef mesh_t::edge_ptr_t edge_ptr_t;
typedef mesh_t::face_ptr_t face_ptr_t;

}


#endif
