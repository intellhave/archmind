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

#ifndef OPENCL_SOLVER_H
#define OPENCL_SOLVER_H

#include "SphereGeometry.h"
#include "OpenCL.h"
#include "Utils.h"

#include <vector>

#ifndef USE_DOUBLE
#define USE_DOUBLE 0
#endif

#if USE_DOUBLE
typedef double scalar_t;
typedef cl_double cl_scalar_t;
typedef cl_double4 cl_scalar4_t;
#else
typedef float scalar_t;
typedef cl_float cl_scalar_t;
typedef cl_float4 cl_scalar4_t;
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
		max_iters(100), 
		max_sp_iters(10000),
		weights(0),
		cpu(false),
		centroid_proj(false)
		{}

	scalar_t target_residual;
	scalar_t spdelta;
	std::size_t workgroup;
	std::size_t max_iters;
	std::size_t max_sp_iters;
	unsigned weights;	//0 - tutte, 1 - conformal
	bool cpu;
	bool centroid_proj;
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

	int m_Times;
	
	mesh_t &m_Input;
	mesh_t &m_Output;

	cl_context m_hContext;
	cl_program m_Program;
	cl_kernel m_Kernel;

	cl_kernel m_KernelSPConvergence;
	cl_kernel m_KernelSRes;

	cl_kernel m_KernelNormalize;
	cl_kernel m_KernelConvergence;
	cl_kernel m_KernelRes;

	cl_command_queue m_CmdQueue;
	std::vector<cl_mem> m_Memobjs;

	cl_scalar4_t *m_Vertices;
	cl_scalar4_t *m_SVertices;
	cl_scalar_t *m_Weights;
	cl_int *m_NVertices;
	std::vector<cl_int> m_Bounds;

	std::size_t m_LocalSize;
	std::size_t m_GlobalSize;
    std::size_t m_ReductionThreads;
    std::size_t m_ReductionBlocks;
	cl_int m_StrideSize;

	Options m_Options;
};

}

#endif