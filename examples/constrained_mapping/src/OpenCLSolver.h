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

#ifndef OPENCL_SOLVER_H
#define OPENCL_SOLVER_H

#define PROFILE_OPENCL false

#include "config.h"
#include <vector>
#include "GeometryTraits.h"

#define __CL_ENABLE_EXCEPTIONS
#include <CL/cl.hpp>


#if CL_USE_DOUBLE
typedef double scalar_t;
typedef cl_double cl_scalar_t;
typedef cl_double2 cl_scalar2_t;
typedef cl_double4 cl_scalar4_t;
typedef cl_double8 cl_scalar8_t;
#else
typedef float scalar_t;
typedef cl_float cl_scalar_t;
typedef cl_float2 cl_scalar2_t;
typedef cl_float4 cl_scalar4_t;
typedef cl_float8 cl_scalar8_t;
#endif

namespace spheremap
{

struct Stats
{
	unsigned elapsed_ms;
	unsigned iterations;
	scalar_t residual;
	bool success;
};

struct Options
{
	Options() : 
		target_residual(1e-07f),
		spdelta(1e-06f),
		workgroup(512),
        scale_iters(100),
		opt_iters(1000), 
		un_iters(1000),
		weights(0),
		device_type(0xFFFFFFFF),
		projection_type(0),
        theta(1.0),
        profile(PROFILE_OPENCL)
		{}

	scalar_t target_residual;
	scalar_t spdelta;
	std::size_t workgroup;
	std::size_t opt_iters;
    std::size_t scale_iters;
    std::size_t un_iters;
	unsigned weights;	//0 - tutte, 1 - conformal
    unsigned energy;    //0 - mips, 1 - knupp
	int device_type;
	int projection_type;
    bool free_boundaries;
    double theta;
    bool profile;
};

class SolverCL 
{
public:
	SolverCL(
		mesh_t &input_mesh,
		mesh_t &output_mesh,
		Options options  );
	~SolverCL();

	bool solve(Stats &solve_stats);

private:
	void compute_weights();
    void compute_mips_weights(std::vector< cl_scalar_t > &angles, std::vector< scalar_t > &factors);
    void circular_projection(int max_iters);
    void untangle_laplace(int max_iters);

	mesh_t &m_Input;
	mesh_t &m_Output;
    std::vector< mesh_t::face_ptr_t > m_BannedFaces;
	std::vector< cl_scalar_t > m_HostPnts;
    cl_scalar_t m_Scale;
    cl_scalar_t m_Diagonal;
    mesh_t::point_t m_Origin;
    
	Options m_Options;
	cl::CommandQueue m_Queue;
	cl::Context m_Context;
	cl::Buffer m_Pnts;
	cl::Buffer m_PntsAux;
	cl::Buffer m_Df;
	cl::Buffer m_Tris;
	cl::Buffer m_Quads;
    cl::Buffer m_Quads_ia;
    cl::Buffer m_Quads_ja;
    cl::Buffer m_Tris_ia;
    cl::Buffer m_Tris_ja;
    cl::Buffer m_TrisAngles;
    cl::Buffer m_TrisFactors;
	std::size_t m_verts_size;
	std::size_t m_tris_size;
	std::size_t m_quads_size;
	std::size_t m_fixed_pnts;
};

}

#endif
