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
#include "SphereGeometry.h"
#include "Io/Io.h"
#include "CommandLine.h"
#include "OpenCLSolver.h"

using namespace spheremap;

int main(int argc, char *argv[])
{
    std::cout << 
        "Parallel Computation of Spherical Parameterizations for Mesh Analysis\n" <<
        "Copyright (C) 2011 Athanasiadis Theodoros and Fudos Ioannis\n\n";

    cmdline_map_t cmd_args;

    if( !parseCommandLine(argc,argv,cmd_args) ) 
        return 1;

    mesh_t input_mesh;
    mesh_t output_mesh;

    if( !arch::io::load_from_file(cmd_args["source"].as<std::string>(), input_mesh) )
    {
        std::cerr << "Failed to open : " << cmd_args["source"].as<std::string>() << "\n";
        return 1;
    }

    Options options;

    options.workgroup = cmd_args["workgroup"].as<int>();
    options.cpu = cmd_args["cpu"].as<int>();
    options.max_iters = cmd_args["max_iters"].as<int>();
    options.spdelta = cmd_args["spdelta"].as<float>();
    options.target_residual = cmd_args["res"].as<float>();
    options.weights = cmd_args["weights"].as<std::string>() == "conformal" ? 1 : 0;
    options.centroid_proj = cmd_args["centroid"].as<int>();

    Stats stats;

    SolverCL solver(input_mesh, output_mesh, options);

    solver.solve(stats);

    std::cout << "Stats\n";
    std::cout << "Time : " << stats.elapsed_ms << " ms\n";
    std::cout << "Iters : " << stats.iterations << "\n";
    std::cout << "Residual : " << stats.residual << "\n";

    arch::io::save_to_file(cmd_args["target"].as<std::string>(), output_mesh);

    return 0;
}