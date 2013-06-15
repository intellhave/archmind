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

#include "PyTypes.h"
#include "PyGeometry.h"
#include "../Geometry/Geometry.h"
#include "../Geometry/Algorithms.h"

#include <boost/python.hpp>
#include <boost/python/stl_iterator.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

using arch::python::real_t;

using namespace arch::geometry;

typedef mesh< default_traits<real_t> > mesh_t;

typedef mesh_t::face_t face_t;
typedef mesh_t::edge_t edge_t;
typedef mesh_t::vertex_t vertex_t;

typedef mesh_t::face_ptr_t face_ptr_t;
typedef mesh_t::edge_ptr_t edge_ptr_t;
typedef mesh_t::vertex_ptr_t vertex_ptr_t;

typedef boost::shared_ptr< mesh_t > mesh_ptr_t;

face_ptr_t triangle( vertex_ptr_t v0, vertex_ptr_t v1, vertex_ptr_t v2)
{
    std::vector< vertex_ptr_t > verts;
    verts.push_back( v0 );
    verts.push_back( v1 );
    verts.push_back( v2 );

    return face_ptr_t(new face_t( verts.begin(), verts.end() ));
}

face_ptr_t quad( vertex_ptr_t v0, vertex_ptr_t v1, vertex_ptr_t v2, vertex_ptr_t v3)
{
    std::vector< vertex_ptr_t > verts;
    verts.push_back( v0 );
    verts.push_back( v1 );
    verts.push_back( v2 );
    verts.push_back( v3 );

    return face_ptr_t(new face_t( verts.begin(), verts.end() ));
}

mesh_ptr_t createmesh()
{
    return mesh_ptr_t( new mesh_t() );
}

template<typename T>
bool py_is_convex(const boost::python::object &o)
{
    typedef boost::python::stl_input_iterator<T> iterator_t;
    iterator_t begin(o), end;
    
    return is_convex( begin, end );
}

//Convert a python list to a c++ iterator range
template<typename T,typename R>
R py_centroid(const boost::python::object &o)
{
    typedef boost::python::stl_input_iterator<T> iterator_t;
    iterator_t begin(o), end;

    return centroid(begin,end);
}

template<typename T, typename R>
R py_area(const boost::python::object &o)
{
    typedef boost::python::stl_input_iterator<T> iterator_t;
    iterator_t begin(o), end;

    return area(begin,end);
}

template<typename T>
T py_circum_center(const boost::python::object &o)
{
    typedef boost::python::stl_input_iterator<T> iterator_t;
    iterator_t begin(o), end;

    return circum_center(begin,end);
}

template<typename T, typename R>
R py_circum_radius(const boost::python::object &o)
{
    typedef boost::python::stl_input_iterator<T> iterator_t;
    iterator_t begin(o), end;

    return circum_radius(begin,end);
}


//Convert a python list to a c++ iterator range
template<typename T>
face_ptr_t py_poly(boost::python::object o)
{
    typedef boost::python::stl_input_iterator<T> iterator_t;

    iterator_t begin(o), end;

    return face_ptr_t( new face_t( begin, end ) );
}

std::vector< vertex_ptr_t > triangulate_ret(const face_ptr_t &f)
{
    std::vector< vertex_ptr_t > tris;
    triangulate( f, tris );

    return tris;
}

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(split_edge_overloads, split_edge, 1, 3)

void arch::python::export_geometry()
{
    using namespace math;
    using namespace boost::python;

    //map the Geometry to a sub-module
    object geometry_module( handle<>(borrowed(PyImport_AddModule("archmind.geometry") ) ) );

    scope().attr("geometry");
    scope geometry_scope = geometry_module;

    def("triangle", &triangle);
    def("quad", &quad);
    def("poly", &py_poly<mesh_t::vertex_ptr_t>);
    def("mesh", &createmesh);

    def("is_free", &is_free< edge_ptr_t > );
    def("is_tjoin", &is_tjoin< edge_ptr_t > );
    
    def("is_valid", &is_valid< face_ptr_t > );
    def("is_valid", &is_valid< edge_ptr_t > );
    def("is_valid", &is_valid< vertex_ptr_t > );

    void (*triangulate_fn)(const face_ptr_t &, std::vector< vertex_ptr_t > &) = 
                &geometry::triangulate< face_ptr_t, std::vector< vertex_ptr_t > >;
    
    def("triangulate", triangulate_fn,
            "Triangulates a simple planar face and save the result in the list", args("f","list"));

    def("triangulate", &triangulate_ret,
            "Triangulates a simple planar face and returns a list of vertices", args("f"));

    def("centroid", &py_centroid< mesh_t::vertex_ptr_t, mesh_t::point_t >);
    def("is_convex", &py_is_convex<  mesh_t::point_t >);
    
    def("circum_center", &py_circum_center< mesh_t::point_t >,
            "Circumcenter of oriented points", args("points"));

    def("circum_radius", &py_circum_radius< mesh_t::point_t, mesh_t::real_t >,
            "Circumradius of oriented points", args("points"));
    
    //Export vector3
    typedef mesh_t::point_t point_t;

    def("area", &py_area< mesh_t::point_t, mesh_t::real_t >, 
            "Signed area of simple polygon", args("points"));
    
    class_< point_t >("vec3")
        .def(init<real_t>())
        .def(init<real_t,real_t,real_t>())
        .def(init<float,float,float>())
        .def(init<const vec3<real_t> &>())
        .def_readwrite("x", &point_t::x)
        .def_readwrite("y", &point_t::y)
        .def_readwrite("z", &point_t::z)
        .def(self + self)
        .def(self - self)
        .def(self += self)
        .def(self -= self)
        .def(self / real_t())
        .def(self * real_t())
        .def(self /= real_t())
        .def(self *= real_t())
        .def( -self )
        //.def(str(self))     // __str__, this should be the correct definition according to the manual
        .def(self_ns::str(self_ns::self))       //but only this works
    ;

    //Export vector3 functions
    def("distance", &math::distance< point_t >);
    def("normalize", &math::normalize< point_t >);

    real_t (*dot_fn)(const point_t &, const point_t &) = &math::dot< real_t >;
    real_t (*adot_fn)(const point_t &, const point_t &) = &math::adot< real_t >;
    point_t (*cross_fn)(const point_t &, const point_t &) = &math::cross< real_t >;
    point_t (*normal_fn)(const point_t &, const point_t &, const point_t &) = &math::normal< real_t >;
    real_t (*mag_fn)(const point_t &) = &math::magnitude< real_t >;
    
    def("dot", dot_fn);
    def("adot", adot_fn);
    def("cross", cross_fn);
    def("normal", normal_fn);
    def("magnitude", mag_fn);

    real_t (*clip_line_to_plane_fn)(const point_t &, const point_t &, const point_t &, const real_t &) = &math::clip_line_to_plane< real_t >;
    bool (*point_in_triangle_fn)(const point_t &, const point_t &, const point_t &, const point_t &) = &math::point_in_triangle< real_t >;

    def("clip_line_to_plane", clip_line_to_plane_fn );
    def("point_in_triangle", point_in_triangle_fn );

    class_< vertex_t, vertex_ptr_t >("vertex", "Mesh vertex")
        .def(init<const point_t &>())
        .def(init<real_t,real_t,real_t>())
        .add_property("faces", range(&vertex_t::faces_begin, &vertex_t::faces_end))
        .add_property("verts", range(&vertex_t::verts_begin, &vertex_t::verts_end))
        .add_property("edges", range(&vertex_t::edges_begin, &vertex_t::edges_end))
        .add_property("verts_size", &vertex_t::verts_size, "Number of vertices")
        .add_property("edges_size", &vertex_t::edges_size, "Number of incident edges")
        .add_property("point", &vertex_t::point, "Vertex coordinates")
        .def(self == self)  //== operator, it is needed by python to do proper comparison
        .def(self != self)  //!= operator
		.def(self < self)
        .def("__hash__", &vertex_t::unique_id)
        .add_property("id", &vertex_t::id, "Index in the mesh")
        .add_property("uid", &vertex_t::unique_id, "Entity unique ID")
    ;

    class_< std::vector< vertex_ptr_t > >("vertex_vec", init<>())
        .def(vector_indexing_suite<std::vector< vertex_ptr_t >, true>())
        .def("__len__", &std::vector< vertex_ptr_t >::size)
        .def("__iter__", iterator<std::vector< vertex_ptr_t > >())
    ;

    class_< edge_t, edge_ptr_t >("edge", "Mesh edge")
        .add_property("verts", range(&edge_t::verts_begin, &edge_t::verts_end), "Adjacent vertices")
        .add_property("verts_size", &edge_t::verts_size, "Number of adjacent vertices" )
        .add_property("faces", range(&edge_t::faces_begin, &edge_t::faces_end), "Incident faces")
        .add_property("faces_size", &edge_t::faces_size, "Number of incident faces")
        .add_property("v0", &edge_t::v0, "First adjacent vertex")
        .add_property("v1", &edge_t::v1, "Second adjacent vertex")
        .def(self == self)
        .def(self != self)
		.def(self < self)
        .def("__hash__", &edge_t::unique_id)
        .add_property("id", &edge_t::id, "Index in the mesh")
        .add_property("uid", &edge_t::unique_id, "Entity unique ID")
    ;

    class_< face_t, face_ptr_t >("face", "Planar mesh face")
        .add_property("points", range(&face_t::points_begin, &face_t::points_end), "Face coordinates")
        .add_property("points_size", &face_t::points_size, "Number of face coordinates")
        .add_property("verts", range(&face_t::verts_begin, &face_t::verts_end),"Adjacent vertices")
        .add_property("verts_size", &face_t::verts_size, "Number of adjacent vertices")
        .add_property("edges", range(&face_t::edges_begin, &face_t::edges_end),"Adjacent edges")
        .add_property("edges_size", &face_t::edges_size,"Number of adjacent edges")
        .add_property("faces", range(&face_t::faces_begin, &face_t::faces_end),"Adjacent faces")
        .def("edge_cw", &face_t::edge_cw, "Returns true if the given edge is visited clockwise by the face", args("e"))
        .def("edge_ccw", &face_t::edge_ccw, "Returns true if the given edge is visited counter-clockwise by the face", args("e"))
        .add_property("normal", &face_t::normal, "Face normal")
        .def(self == self)
        .def(self != self)
		.def(self < self)
        .def("__hash__", &face_t::unique_id)
        .add_property("id", &face_t::id, "Index in the mesh")
        .add_property("uid", &face_t::unique_id, "Entity unique ID")
    ;
    
    class_< mesh_t, boost::shared_ptr< mesh_t > >("mesh","Non-manifold surface mesh")
        .def("set_point", &mesh_t::set_point, "Changes the coordinates of a vertex", args("v","point"))
        .def("add_vertex", &mesh_t::add_vertex, "Adds a free vertex to the mesh and returns the index", args("v") )
        .def("remove_vertex", &mesh_t::remove_vertex)
        .def("add_face", &mesh_t::add_face)
        .def("remove_face", &mesh_t::remove_face)
        .def("split_edge", &mesh_t::split_edge, split_edge_overloads())
        .def("join_edge", &mesh_t::join_edge)
        .def("join_face", &mesh_t::join_face) 
        .def("split_face", &mesh_t::split_face)
		.def("flip_face", &mesh_t::flip_face, "Flips the orientation of a face", args("f"))
        .add_property("verts", range(&mesh_t::verts_begin, &mesh_t::verts_end))
        .add_property("verts_size", &mesh_t::verts_size)
        .add_property("edges", range(&mesh_t::edges_begin, &mesh_t::edges_end))
        .add_property("edges_size", &mesh_t::edges_size)
        .add_property("faces", range(&mesh_t::faces_begin, &mesh_t::faces_end))
        .add_property("faces_size", &mesh_t::faces_size)
        ;   
    
}
