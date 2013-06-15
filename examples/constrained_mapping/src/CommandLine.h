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

#ifndef CMDLINE_H
#define CMDLINE_H

#include <boost/program_options.hpp>

namespace parameterization
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
      ("map", value<string>()->default_value(""), "source pinned vertices file")
      ("target", value<string>(), "target model file")
      ("opt_iters", value<std::size_t>()->default_value(1000), "non linear optimizer max iterations")
      ("un_iters", value<std::size_t>()->default_value(1000), "Untangling process max iterations")
      ("scale_iters", value<std::size_t>()->default_value(300), "scale iterations")
      ("device", value<int>()->default_value(0xFFFFFFFF), "set to 2 to run the solver on the cpu or 4 to run on the gpu")
      ("res", value<float>()->default_value(1e-07f), "target residual")
      ("workgroup", value<int>()->default_value(512), "OpenCL workgroup size")
      ("proj", value<int>()->default_value(2), "Initial projection (0 = planar, 1 = circular, 2 = use uv)")
      ("free", value<int>()->default_value(1), "Free boundaries")
      ("type", value<string>()->default_value("isometric"), "type of parameterization (mips,isometric,smooth)")
      ("ps", value<int>()->default_value(0), "Export post script (1 = gnuplot data, 2 = plain mesh, 3 = mesh with area deformation)")
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

    if( vm.count("source") || vm.count("target")  )
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
