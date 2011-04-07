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

#ifndef CMDLINE_H
#define CMDLINE_H

#include <boost/program_options.hpp>

namespace spheremap
{

typedef boost::program_options::variables_map cmdline_map_t;

bool parseCommandLine(int argc, char **argv, cmdline_map_t &vm)
{
	using namespace boost::program_options;
	using namespace std;

	//Set the available options
	options_description desc("Allowed Command Line Arguments");
		desc.add_options()
			("help", "produce help message")
			("source", value<string>(), "source model file")
			("target", value<string>(), "target model file")
			("max_iters", value<int>()->default_value(100), "max iterations")
			("weights", value<string>()->default_value("tutte"), "type of weights (tutte,conformal)")
			("cpu", value<int>()->default_value(0), "set to 1 to run the solver on the cpu")
			("res", value<float>()->default_value(1e-07f), "target residual")
			("spdelta", value<float>()->default_value(1e-06f), "saddle point convergence delta")
			("workgroup", value<int>()->default_value(512), "OpenCL workgroup size")
			("centroid", value<int>()->default_value(0), "set to 1 to use a centroidal projection for the initial solution")
			;

	try
	{
		store(parse_command_line(argc,argv,desc), vm);
		notify(vm);

		if(vm.count("help"))
		{
			cout << desc << endl;		//output help
			return false;
		}
	}
	catch(exception &e)
	{
		cerr << e.what() << endl;
		return false;
	}

	if( vm.count("source") && vm.count("target") )
	{
		return true;
	}
	else
	{
		std::cout << "Type : " << argv[0] << " --help for a list of options\n";
		return false;
	}
}

}

#endif