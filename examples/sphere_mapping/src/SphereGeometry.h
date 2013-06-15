/*
  Parallel Computation of Spherical Parameterizations for Mesh Analysis
  Copyright (C) 2011 Athanasiadis Theodoros and Fudos Ioannis

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
struct sphere_vertex : arch::geometry::vertex<Traits>
{
    typedef typename Traits::point_t point_t;
    typedef typename Traits::real_t real_t;

	sphere_vertex(const point_t &co) : arch::geometry::vertex<Traits>(co) {}
	sphere_vertex(const real_t &x, const real_t &y, const real_t &z) : arch::geometry::vertex<Traits>(x,y,z) {}

	std::vector< float > Weights;
};

struct sphere_traits : 
	arch::geometry::traits<
        float,
        sphere_vertex< sphere_traits >,
		arch::geometry::edge< sphere_traits >,
        arch::geometry::face< sphere_traits >,
        arch::geometry::mesh< sphere_traits > 
    >
{
};

typedef arch::geometry::mesh<sphere_traits> mesh_t;
typedef mesh_t::point_t point_t;

typedef mesh_t::vertex_t vertex_t;
typedef mesh_t::edge_t edge_t;
typedef mesh_t::face_t face_t;
typedef mesh_t::vertex_ptr_t vertex_ptr_t;
typedef mesh_t::edge_ptr_t edge_ptr_t;
typedef mesh_t::face_ptr_t face_ptr_t;

}


#endif
