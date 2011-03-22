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

#include "OpenCLSolver.h"
#include <iostream>
#include <vector>

#if defined (_WIN32)	
	#include <windows.h>
	#include <mmsystem.h>
	#pragma comment(lib, "Winmm.lib")
#else
	#include <sys/time.h>

	unsigned int get_ticks()
	{
		struct timeval t;

		gettimeofday(&t, 0);
		return (t.tv_sec*1000) + (t.tv_usec/1000); //bring down number to miliseconds
	}

	#define timeGetTime() (get_ticks())
#endif

using namespace std;
using namespace arch::math;
using namespace arch::geometry;
using namespace spheremap;

spheremap::SolverCL::SolverCL( 
	mesh_t &input_mesh,
	mesh_t &ouput_mesh,
	Options options) : m_Input(input_mesh), m_Output(ouput_mesh), m_Options(options)
{
	std::size_t n = input_mesh.vertices_size();
	std::size_t m = input_mesh.faces_size();

	m_LocalSize = 1024;
	m_GlobalSize = cl::roundUp(m_LocalSize,n);

	//Vertices
	m_Vertices = new cl_float4[n];
	m_SVertices = new cl_float4[n];		//Output

	point_t c = options.centroid_proj ? centroid( input_mesh.vertices() ) : point_t(0.0);

	//Copy the vertices
	std::size_t i = 0;
	for( mesh_t::vertex_iterator_t v = input_mesh.vertices_begin(); 
		v != input_mesh.vertices_end(); ++v, ++i )
	{
		//const point_t &p = (*v)->point();
		point_t p = normalize( (*v)->point() - c );

		m_Vertices[i].s[0] = p.x;
		m_Vertices[i].s[1] = p.y;
		m_Vertices[i].s[2] = p.z;
		m_Vertices[i].s[3] = 0.0f;
	}

	//Copy the neighboring information
	for( std::size_t i = 0; i < n; i++ )
	{
		m_Bounds.push_back( i == 0 ? 0 : m_Bounds.back() );
		m_Bounds.push_back( m_Bounds.back() + input_mesh.vertices()[i]->edges_size() );
	}

	//std::cout << "Computing conformal weights...\n";
	compute_weights();

	m_NVertices = new cl_uint[m_Bounds.back()];
	m_Weights = new cl_float[m_Bounds.back()];

	cl_uint *nverts = m_NVertices;
	cl_float *nweights = m_Weights;

	//Copy the neighboring information
	for( mesh_t::vertex_iterator_t v = input_mesh.vertices_begin(); 
		v != input_mesh.vertices_end(); ++v )
	{
		std::size_t wi = 0;
		for( mesh_t::vertex_t::vertex_iterator_t vv = (*v)->vertices_begin(); 
			vv != (*v)->vertices_end(); ++vv)
		{
			*nverts++ = (*vv)->get_id();
			*nweights++ = (*v)->Weights[wi++];
		}
	}
	
	cl_int err = 0;

	cl_platform_id platform = NULL;
	std::vector< cl_platform_id > platforms;

	cl::get_platforms(platforms);

	for( std::size_t i = 0; i < platforms.size(); ++i )
	{
		std::string name;
		cl::platform_info(platforms[i],CL_PLATFORM_NAME,name);
		if( (m_Options.cpu && name.find("ATI") != string::npos) || (!m_Options.cpu && name.find("NVIDIA") != string::npos)  )
		{
			std::cout << "Using : " << name << '\n';
			platform = platforms[i];
			break;
		}
	}

	if( platform == NULL )
	{	
		std::cerr << "Error : Failed to find platform\n"; 
		return;
	}

	cl_device_id device;
	
	if( m_Options.cpu )
		clGetDeviceIDs(platform, CL_DEVICE_TYPE_CPU, 1, &device, 0);
	else
		clGetDeviceIDs(platform, CL_DEVICE_TYPE_GPU, 1, &device, 0);

	// create the OpenCL context on a GPU device
	m_hContext = clCreateContext( 0, 1, &device, 0, 0, &err );

	if( err != CL_SUCCESS )
	{
		std::cerr << "Error in clCreateContextFromType : " << err << std::endl;
		return;
	}

	// get the list of GPU devices associated with context
	size_t nContextDescriptorSize;
	clGetContextInfo(m_hContext, CL_CONTEXT_DEVICES, 0, 0, &nContextDescriptorSize);
	cl_device_id *aDevices = (cl_device_id *)malloc(nContextDescriptorSize);
	clGetContextInfo(m_hContext, CL_CONTEXT_DEVICES, nContextDescriptorSize, aDevices, 0);

	clGetDeviceInfo(aDevices[0], CL_DEVICE_MAX_WORK_GROUP_SIZE, sizeof(std::size_t), &m_LocalSize, 0);
	m_LocalSize = m_Options.workgroup;
	std::cout << "OpenCL : work group size : " << m_LocalSize << '\n';

	m_GlobalSize = cl::roundUp(m_LocalSize,n);

	// create a command-queue
	m_CmdQueue = clCreateCommandQueue(m_hContext, aDevices[0], 0, 0);
	
	// allocate the buffer memory objects
	//Vertices
	m_Memobjs.push_back(clCreateBuffer(m_hContext, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(cl_float4)*n, m_Vertices, NULL));
	//Constraints
	m_Memobjs.push_back(clCreateBuffer(m_hContext, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(cl_float4)*n, m_Vertices, NULL));
	//Bounds
	m_Memobjs.push_back(clCreateBuffer(m_hContext, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(cl_uint)*m_Bounds.size(), &m_Bounds[0], NULL));
	//Indices
	m_Memobjs.push_back(clCreateBuffer(m_hContext, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(cl_uint)*m_Bounds.back(), m_NVertices, NULL));
	//Weights
	m_Memobjs.push_back(clCreateBuffer(m_hContext, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(cl_float)*m_Bounds.back(), m_Weights, NULL));
	m_Memobjs.push_back(clCreateBuffer(m_hContext, CL_MEM_READ_WRITE, sizeof(cl_float4)*n, NULL, NULL));
	m_Memobjs.push_back(clCreateBuffer(m_hContext, CL_MEM_WRITE_ONLY, sizeof(cl_float)*n, NULL, NULL));
	m_Memobjs.push_back(clCreateBuffer(m_hContext, CL_MEM_WRITE_ONLY, sizeof(cl_float), NULL, NULL));

	std::string program_source;
	if( loadTextFile("CL/nl_laplacesolver.cl",program_source) )
		std::cout << "OpenCL file" << " successfully loaded!" << std::endl;

	const char *source = program_source.c_str();

	// create the program
	m_Program = clCreateProgramWithSource(m_hContext, 1, &source, 0, &err);

	if( err != CL_SUCCESS )
	{
		std::cerr << "Error in clCreateProgramWithSource!" << std::endl;
		return;
	}

	// build the program
	err = clBuildProgram(m_Program, 0, NULL, NULL, NULL, NULL);

	if( err != CL_SUCCESS )
	{
		std::cerr << "Error in clBuildProgram!" << std::endl;
		std::cerr << cl::buildlog(m_Program,aDevices[0]) << std::endl;
		return;
	}

	if( m_Options.weights == 0 )
	{
		m_Kernel = clCreateKernel(m_Program, "solve_equal", NULL);
		m_KernelSPConvergence = clCreateKernel(m_Program, "solve_equal_res", NULL);
		m_KernelConvergence = clCreateKernel(m_Program, "convergence_equal_res", NULL);
	}
	else
	{
		m_Kernel = clCreateKernel(m_Program, "solve_conformal", NULL);
		m_KernelSPConvergence = clCreateKernel(m_Program, "solve_conformal_res", NULL);
		m_KernelConvergence = clCreateKernel(m_Program, "convergence_conformal_res", NULL);
	}

	m_KernelSRes = clCreateKernel(m_Program, "lmax_residual", NULL);
	m_KernelNormalize = clCreateKernel(m_Program, "normalize_solution", NULL);
	m_KernelRes = clCreateKernel(m_Program, "l2_residual", NULL);
}

spheremap::SolverCL::~SolverCL()
{
	//clean up
	clReleaseKernel( m_Kernel );
	clReleaseKernel( m_KernelSPConvergence );
	clReleaseKernel( m_KernelSRes );
	clReleaseKernel( m_KernelNormalize );
	clReleaseKernel( m_KernelRes );
	clReleaseKernel( m_KernelConvergence );

	clReleaseProgram( m_Program );
	clReleaseCommandQueue( m_CmdQueue );
	clReleaseContext( m_hContext );
	
	for( std::size_t i = 0; i < m_Memobjs.size(); ++i )
		clReleaseMemObject( m_Memobjs[i] );

	delete [] m_Vertices;
	delete [] m_SVertices;
	delete [] m_NVertices;
	delete [] m_Weights;
}

bool spheremap::SolverCL::solve(spheremap::Stats &solve_stats)
{
	using std::max;

	solve_stats.elapsed_ms = timeGetTime();
	solve_stats.success = true;
	solve_stats.iterations = 0;

	std::size_t n = m_Input.vertices_size();
	std::size_t m = m_Input.faces_size();

	cl_mem src,dst,weights,constraints,bounds,indices,tmp_res;

	src = m_Memobjs[0];		//Vertices
	constraints = m_Memobjs[1];		//Constraints
	bounds = m_Memobjs[2];		//Bounds
	indices = m_Memobjs[3];		//Indices
	weights = m_Memobjs[4];		//Weights
	dst = m_Memobjs[5];
	tmp_res = m_Memobjs[6];

	//set kernel arguments
	cl_int err = clSetKernelArg(m_Kernel, 0, sizeof(cl_mem), (void *)&src );
	err |= clSetKernelArg(m_Kernel, 1, sizeof(cl_mem), (void *)&constraints );
	err |= clSetKernelArg(m_Kernel, 2, sizeof(cl_mem), (void *)&bounds );
	err |= clSetKernelArg(m_Kernel, 3, sizeof(cl_mem), (void *)&indices );
	err |= clSetKernelArg(m_Kernel, 4, sizeof(cl_mem), (void *)&dst );
	err |= clSetKernelArg(m_Kernel, 5, sizeof(size_t), &n );

	if( m_Options.weights == 1 )	//conformal weights
		err |= clSetKernelArg(m_Kernel, 6, sizeof(cl_mem), (void *)&weights );

	err |= clSetKernelArg(m_KernelSPConvergence, 0, sizeof(cl_mem), (void *)&src );
	err |= clSetKernelArg(m_KernelSPConvergence, 1, sizeof(cl_mem), (void *)&constraints );
	err |= clSetKernelArg(m_KernelSPConvergence, 2, sizeof(cl_mem), (void *)&bounds );
	err |= clSetKernelArg(m_KernelSPConvergence, 3, sizeof(cl_mem), (void *)&indices );
	err |= clSetKernelArg(m_KernelSPConvergence, 4, sizeof(cl_mem), (void *)&dst );
	err |= clSetKernelArg(m_KernelSPConvergence, 5, sizeof(size_t), &n );
	err |= clSetKernelArg(m_KernelSPConvergence, 6, sizeof(cl_mem), (void *)&tmp_res );

	if( m_Options.weights == 1 )	//conformal weights
		err |= clSetKernelArg(m_KernelSPConvergence, 7, sizeof(cl_mem), (void *)&weights );
	
	err |= clSetKernelArg(m_KernelConvergence, 1, sizeof(cl_mem), (void *)&bounds );
	err |= clSetKernelArg(m_KernelConvergence, 2, sizeof(cl_mem), (void *)&indices );
	err |= clSetKernelArg(m_KernelConvergence, 3, sizeof(size_t), &n );
	err |= clSetKernelArg(m_KernelConvergence, 4, sizeof(cl_mem), (void *)&tmp_res );

	if(  m_Options.weights == 1 )
		err |= clSetKernelArg(m_KernelConvergence, 5, sizeof(cl_mem), (void *)&weights );

	if( err != CL_SUCCESS )
	{
		std::cerr << "Some error occured during the setting of arguments : " << err << std::endl;
		return false;
	}

	cl_float last_res = 1000000.0f;
	cl_float res, linear_res;
	std::size_t one = 1;

	for( int i = 0; i < m_Options.max_iters; ++i )
	{
		solve_stats.iterations++;

		for( int j = 1; j < m_Options.max_sp_iters; ++j )
		{
			// execute the normal kernel
			if( j % 1000 != 0 )
			{
				// execute kernel 
				err = clEnqueueNDRangeKernel(m_CmdQueue, m_Kernel, 1, NULL, &m_GlobalSize, &m_LocalSize, 0, NULL, NULL);

				if( err != CL_SUCCESS )
				{
					std::cerr << "Some error occured during the execution of the kernel : " << err << std::endl;
					return false;
				}

				//swap buffers
				std::swap(dst,src);

				//reset arguments
				clSetKernelArg(m_Kernel, 0, sizeof(cl_mem), (void *)&src );
				clSetKernelArg(m_Kernel, 4, sizeof(cl_mem), (void *)&dst );
			}
			// execute the kernel that also calculates the residual
			else
			{
				// execute kernel 
				err = clEnqueueNDRangeKernel(m_CmdQueue, m_KernelSPConvergence, 1, NULL, &m_GlobalSize, &m_LocalSize, 0, NULL, NULL);

				//swap buffers
				std::swap(dst,src);

				//reset arguments
				clSetKernelArg(m_Kernel, 0, sizeof(cl_mem), (void *)&src );
				clSetKernelArg(m_Kernel, 4, sizeof(cl_mem), (void *)&dst );

				//Finally compute the residual (Lmax)
				clSetKernelArg(m_KernelSRes, 0, sizeof(size_t), &n );
				clSetKernelArg(m_KernelSRes, 1, sizeof(cl_mem), &tmp_res );
				clSetKernelArg(m_KernelSRes, 2, sizeof(cl_mem), &m_Memobjs.back());

				err |= clEnqueueNDRangeKernel(m_CmdQueue, m_KernelSRes, 1, NULL, &one, &one, 0, NULL, NULL);

				if( err != CL_SUCCESS )
				{
					std::cerr << "Some error occured during the execution of the kernel : " << err << std::endl;
					return false;
				}

				//Read back the residual value
				err |= clEnqueueReadBuffer(m_CmdQueue, m_Memobjs.back(), CL_TRUE, 0, sizeof(cl_float), &linear_res, 0, NULL, NULL);

				if( linear_res < m_Options.spdelta )
					break;				
			}
		}

		//Normalize the solution and update the constraints
		clSetKernelArg(m_KernelNormalize, 0, sizeof(cl_mem), (void *)&src);			//vertices
		clSetKernelArg(m_KernelNormalize, 1, sizeof(size_t), (void *)&n);			//vertices
		clSetKernelArg(m_KernelNormalize, 2, sizeof(cl_mem), (void *)&constraints);	//constraints
		clSetKernelArg(m_KernelNormalize, 3, sizeof(cl_mem), (void *)&dst);			//vertices

		err = clEnqueueNDRangeKernel(m_CmdQueue, m_KernelNormalize, 1, NULL, &m_GlobalSize, &m_LocalSize, 0, NULL, NULL);

		//swap buffers
		std::swap(dst,src);

		//Calculate the residual for each vertex
		clSetKernelArg(m_KernelConvergence, 0, sizeof(cl_mem), (void *)&src );
		
		err = clEnqueueNDRangeKernel(m_CmdQueue, m_KernelConvergence, 1, NULL, &m_GlobalSize, &m_LocalSize, 0, NULL, NULL);

		//Finally compute the residual (L2)
		clSetKernelArg(m_KernelRes, 0, sizeof(size_t), &n );
		clSetKernelArg(m_KernelRes, 1, sizeof(cl_mem), &tmp_res );
		clSetKernelArg(m_KernelRes, 2, sizeof(cl_mem), &m_Memobjs.back());
		
		err = clEnqueueNDRangeKernel(m_CmdQueue, m_KernelRes, 1, NULL, &one, &one, 0, NULL, NULL);

		//Read back the residual value
		err = clEnqueueReadBuffer(m_CmdQueue, m_Memobjs.back(), CL_TRUE, 0, sizeof(cl_float), &res, 0, NULL, NULL);
		
		std::cout << i << " - res : " << res << '\n';

		if( last_res < res || res < m_Options.target_residual )
			break;

		last_res = res;
	}

	//Read the final coordinates
	err = clEnqueueReadBuffer(m_CmdQueue, src, CL_TRUE, 0, n*sizeof(cl_float4), m_SVertices, 0, NULL, NULL);

	solve_stats.elapsed_ms = timeGetTime() - solve_stats.elapsed_ms;

	
	for( std::size_t i = 0; i < n; ++i )
	{
		m_Output.add_vertex( vertex_ptr_t( new vertex_t(
			m_SVertices[i].s[0],m_SVertices[i].s[1],m_SVertices[i].s[2] ) ) );
	}

	for( mesh_t::face_iterator_t f = m_Input.faces_begin(); f != m_Input.faces_end(); ++f )
	{
		std::vector< vertex_ptr_t > verts;
		std::vector< vertex_ptr_t > indices((*f)->vertices_begin(),(*f)->vertices_end());

		for( std::size_t i = 0 ; i < indices.size(); ++i )
			verts.push_back( m_Output.vertices()[ indices[i]->get_id() ] );

		m_Output.add_face( face_ptr_t( new face_t(verts.begin(),verts.end()) ) );
	}

	solve_stats.residual = res;

	return true;
}

void spheremap::SolverCL::compute_weights()
{
	//Compute Conformal or Barycentric weights for all the edges
	for( mesh_t::edge_iterator_t e = m_Input.edges_begin();
		e != m_Input.edges_end(); ++e )
	{
		std::vector< vertex_ptr_t > verts((*e)->vertices_begin(),(*e)->vertices_end());

		//Find the two oposing vertices
		for( edge_t::face_iterator_t ef = (*e)->faces_begin(); 
			ef != (*e)->faces_end(); ++ef )
		{
			for( face_t::vertex_iterator_t efv = (*ef)->vertices_begin();
				efv != (*ef)->vertices_end(); ++efv )
			{
				if( *efv != verts[0] && *efv != verts[1] )
					verts.push_back(*efv);
			}
		}

		float w = 1.0;

		if( m_Options.weights == 1 )	//Conformal weights
		{
			float a0 = acos( dot( 
				normalize( verts[2]->point() - verts[0]->point() ), 
				normalize( verts[2]->point() - verts[1]->point() ) ) );

			float a1 = acos( dot( 
				normalize( verts[3]->point() - verts[0]->point() ), 
				normalize( verts[3]->point() - verts[1]->point() ) ) );

			//Clamp the weights between 5 and 85 degrees
			if( a0 < 0.087266462599716474f ) a0 = 0.087266462599716474f;
			else if( a0 > 1.4835298641951802f ) a0 = 1.4835298641951802f;

			if( a1 < 0.087266462599716474f ) a1 = 0.087266462599716474f;
			else if( a1 > 1.4835298641951802f ) a1 = 1.4835298641951802f;

			w = (1.0 / tan(a0) + 1.0 / tan(a1));
		}

		//add the weight to the two vertices of this edge
		verts[0]->Weights.push_back(w);
		verts[1]->Weights.push_back(w);
	}

	//Normalize the weights of all the vertices
	for( mesh_t::vertex_iterator_t v = m_Input.vertices_begin();
		v != m_Input.vertices_end(); ++v )
	{
		float sum = 0.0;
		for( std::size_t j = 0; j < (*v)->Weights.size(); ++j )
			sum += (*v)->Weights[j];

		for( std::size_t j = 0; j < (*v)->Weights.size(); ++j )
			(*v)->Weights[j] /= sum;
	}
}