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
#include <boost/python.hpp>
#include <iostream>

//#include <cln/float_io.h>

void arch::python::export_types()
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
