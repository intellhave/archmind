#include "PyInterface.h"
#include "PyLoader.h"

#include "../Loaders/Loader.h"
#include "../Geometry/Geometry.h"

void python::export_loader()
{
	using namespace loader;
	using namespace boost::python;

	//map the Geometry to a sub-module
	object geometry_module( handle<>(borrowed(PyImport_AddModule("archmind.io") ) ) );

	scope().attr("io");
	scope geometry_scope = geometry_module;

	def("load_from_file", &load_from_file< geometry::mesh<> >);
	def("save_to_file", &save_to_file< geometry::mesh<> >);
}
