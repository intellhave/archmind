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

#include "PyInterface.h"
#include "PyLoader.h"

#include "../Io/Io.h"
#include "../Geometry/Geometry.h"

void arch::python::export_loader()
{
    using namespace arch::io;
    using namespace boost::python;

    //map the Geometry to a sub-module
    object geometry_module( handle<>(borrowed(PyImport_AddModule("archmind.io") ) ) );

    scope().attr("io");
    scope geometry_scope = geometry_module;

    def("load_from_file", &load_from_file< arch::geometry::mesh<> >);
    def("save_to_file", &save_to_file< arch::geometry::mesh<> >);
}
