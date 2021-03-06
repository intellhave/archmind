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
#include "PyGeometry.h"
#include "PyLoader.h"
#include <Python.h>
#include <locale>

using namespace arch::python;

BOOST_PYTHON_MODULE(archmind)
{
    using namespace boost::python;

    //specify that this module is actually a package
    //object package = scope();
    //package.attr("__path__") = "archmind";
	object archmind_module( handle<>(borrowed(PyImport_AddModule("archmind") ) ) );

    export_types();
    export_geometry();
    export_loader();
}

Interface::Interface(int argc, char **argv) : m_Argc(argc)
{
    using namespace boost::python;

    Py_Initialize();

	//Convert to wchar_t**
	m_Argv = new wchar_t*[argc];
	for( int i = 0; i < argc; ++i )
	{
		m_Argv[i] = new wchar_t[256];
		std::mbstowcs(m_Argv[i],argv[i],255);
	}

    PySys_SetArgv(argc,m_Argv);
	
    object main_module(handle<>(borrowed(PyImport_AddModule("__main__"))));
    m_MainNamespace = main_module.attr("__dict__");
}

Interface::~Interface()
{
    Py_Finalize();

	for( int i = 0; i < m_Argc; ++i )
		delete [] m_Argv[i];

	delete [] m_Argv;
}

std::string Interface::last_error()const
{
    return m_ErrorString;
}

bool Interface::run_script(const std::string &filename)
{
    using namespace boost::python;

    m_ErrorString = "";

    try
    {
		init_module_archmind();
      
        PyRun_SimpleString("import io");
        PyRun_SimpleString("import sys");
        PyRun_SimpleString("sys.stderr = io.StringIO()");

        object result = exec_file(filename.c_str(), m_MainNamespace, m_MainNamespace);
    }
    catch(error_already_set const &)
    {
        PyErr_Print();

        object sys(handle<>(PyImport_ImportModule("sys")));
        
        object err = sys.attr("stderr");
        m_ErrorString = extract<std::string>(err.attr("getvalue")());

        return false;
    }

    return true;
}


