#include "PyTypes.h"
#include <boost/python.hpp>
#include <iostream>

//#include <cln/float_io.h>

void python::export_types()
{
    using namespace boost::python;

	//map the Geometry to a sub-module
	object geometry_module( handle<>(borrowed(PyImport_AddModule("archmind.types") ) ) );

	scope().attr("types");
	scope geometry_scope = geometry_module;

#if 0
    class_< python::real_t > ("real")
        .def(init<const char *>())
        .def(init<real_t>())
        .def(self + self)
        .def(self += self)
        .def(self - self)
        .def(self -= self)
        .def(self * self)
        .def(self *= self)
        .def(self / self)
        .def(self /= self)
        //.def(str(self))
        //.def(float_(self))
        ;
#endif
}
