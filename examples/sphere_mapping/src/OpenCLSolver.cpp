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
#include <boost/foreach.hpp>

#if defined (_WIN32)	
#include <windows.h>
#include <mmsystem.h>
#pragma comment(lib, "Winmm.lib")
#else

#include <sys/time.h>
#include <stdio.h>

unsigned int get_ticks()
{
    struct timeval t;

    gettimeofday(&t, 0);
    return (t.tv_sec*1000) + (t.tv_usec/1000); //bring down number to miliseconds
}

#define timeGetTime() (get_ticks())
#endif

#define foreach BOOST_FOREACH

using namespace std;
using namespace arch::math;
using namespace arch::geometry;
using namespace spheremap;

spheremap::SolverCL::SolverCL( 
                              mesh_t &input_mesh,
                              mesh_t &ouput_mesh,
                              Options options) : m_Input(input_mesh), m_Output(ouput_mesh), m_Options(options)
{
    using std::min;

    std::size_t n = input_mesh.verts_size();

    //Vertices
    m_Vertices = new cl_scalar4_t[n];
    m_SVertices = new cl_scalar4_t[n];		//Output

    point_t c = options.centroid_proj ? centroid( input_mesh.verts() ) : point_t(0.0);

    //Copy the vertices
    std::size_t i = 0;

    foreach( mesh_t::vertex_ptr_t v, input_mesh.verts() )
    {
        point_t p = normalize( v->point() - c );

        m_Vertices[i].s[0] = p.x;
        m_Vertices[i].s[1] = p.y;
        m_Vertices[i].s[2] = p.z;
        m_Vertices[i].s[3] = 0.0f;     
        ++i;
    }

    //Copy the neighboring information
    m_StrideSize = 0;      //distance between successive neighbors 
    for( std::size_t i = 0; i < n; i++ )
    {
        std::size_t wi = input_mesh.verts()[i]->edges_size();
        m_Bounds.push_back( wi );
        if( wi > m_StrideSize ) m_StrideSize = wi;
    }

    compute_weights();

    m_NVertices = new cl_int[ n * m_StrideSize ];
    m_Weights = new cl_scalar_t[ n * m_StrideSize ];

    cl_int *nverts = m_NVertices;
    cl_scalar_t *nweights = m_Weights;

    //Copy the neighboring information
    foreach( mesh_t::vertex_ptr_t v, input_mesh.verts() )
    {
        std::size_t counter = 0;
        std::size_t wi = 0;
        std::size_t last = 0;

        foreach( mesh_t::vertex_t::vertex_ptr_t vv, v->verts() )
        {
            last = counter;
            nverts[ v->get_id() + counter ] = vv->get_id();
            nweights[ v->get_id() + counter ] = v->Weights[wi++];
            counter += n;
        }

        //pad the rest with 0 weights
        for( ; wi < m_StrideSize; ++wi )
        {
            nverts[ v->get_id() + counter ] = nverts[ last ];
            nweights[ v->get_id() + counter ] = 0.0;
            counter += n;
        }	
    }

    cl_int err = 0;

    cl_platform_id platform = NULL;
    std::vector< cl_platform_id > platforms;

    cl::get_platforms(platforms);
    bool valid_platform = false;

    for( std::size_t i = 0; i < platforms.size(); ++i )
    {
        std::string name;
        cl::platform_info(platforms[i],CL_PLATFORM_NAME,name);

        if( (name.find("AMD") != string::npos) || (name.find("NVIDIA") != string::npos)  )
        {
            platform = platforms[i];

            cl_device_id device;
            err = clGetDeviceIDs(platform, m_Options.cpu ? CL_DEVICE_TYPE_CPU : CL_DEVICE_TYPE_GPU, 1, &device, 0);

            if( err != CL_SUCCESS ) continue;

            // create the OpenCL context on a GPU device
            m_hContext = clCreateContext( 0, 1, &device, 0, 0, &err );

            if( err == CL_SUCCESS ) 
            {
                std::cout << "Using : " << name << '\n';
                valid_platform = true;
                break;
            }
        }
    }

    if( !valid_platform )
        throw std::runtime_error("Failed to find platform\n");

    // get the list of GPU devices associated with context
    size_t nContextDescriptorSize;
    clGetContextInfo(m_hContext, CL_CONTEXT_DEVICES, 0, 0, &nContextDescriptorSize);
    cl_device_id *aDevices = (cl_device_id *)malloc(nContextDescriptorSize);
    clGetContextInfo(m_hContext, CL_CONTEXT_DEVICES, nContextDescriptorSize, aDevices, 0);

    clGetDeviceInfo(aDevices[0], CL_DEVICE_MAX_WORK_GROUP_SIZE, sizeof(std::size_t), &m_LocalSize, 0);
    m_LocalSize = min( m_LocalSize, m_Options.workgroup );
    std::cout << "OpenCL : work group size : " << m_LocalSize << '\n';

    m_ReductionThreads = 256;
    m_ReductionBlocks = (n + m_ReductionThreads - 1) / m_ReductionThreads;

    m_GlobalSize = cl::round_up(m_LocalSize,n);

    // create a command-queue
    m_CmdQueue = clCreateCommandQueue(m_hContext, aDevices[0], 0, 0);

    // allocate the buffer memory objects

    //Vertices
    m_Memobjs.push_back(clCreateBuffer(m_hContext, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(cl_scalar4_t) * n, m_Vertices, NULL));   
    //Constraints
    m_Memobjs.push_back(clCreateBuffer(m_hContext, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(cl_scalar4_t) * n, m_Vertices, NULL));
    //Indices
    m_Memobjs.push_back(clCreateBuffer(m_hContext, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(cl_int) * n * m_StrideSize, m_NVertices, NULL));   
    //Weights
    m_Memobjs.push_back(clCreateBuffer(m_hContext, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(cl_scalar_t) * n * m_StrideSize, m_Weights, NULL));  
    //Output
    m_Memobjs.push_back(clCreateBuffer(m_hContext, CL_MEM_READ_WRITE, sizeof(cl_scalar4_t)*n, NULL, NULL));  
    //Temp Residual
    m_Memobjs.push_back(clCreateBuffer(m_hContext, CL_MEM_WRITE_ONLY, sizeof(cl_scalar_t)*n, NULL, NULL)); 
    //Bounds
    m_Memobjs.push_back(clCreateBuffer(m_hContext, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(cl_int)*m_Bounds.size(), &m_Bounds[0], NULL)); 
    //Residual
    m_Memobjs.push_back(clCreateBuffer(m_hContext, CL_MEM_READ_WRITE, sizeof(cl_scalar_t)*m_ReductionBlocks, NULL, NULL));
   
    std::string program_source;
    if( loadTextFile("CL/solver.cl",program_source) )
        std::cout << "OpenCL file" << " loaded" << std::endl;

    const char *source = program_source.c_str();

    // create the program
    m_Program = clCreateProgramWithSource(m_hContext, 1, &source, 0, &err);
    if( err != CL_SUCCESS )
    {
        throw std::runtime_error("Error in clCreateProgramWithSource!");
    }

    // build the program
    char cloptions[255];
    sprintf(cloptions, "-cl-strict-aliasing -D STRIDE_SIZE=%d -D OMEGA=1.0f -D ONE_MINUS_OMEGA=0.0f -D USE_DOUBLE=%d", m_StrideSize, USE_DOUBLE );
    err = clBuildProgram(m_Program, 0, NULL, cloptions, NULL, NULL);

    if( err != CL_SUCCESS )
    {
        std::cerr << cl::buildlog(m_Program,aDevices[0]) << std::endl;
        throw std::runtime_error("Error in clBuildProgram!");
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

    cl_int n = m_Input.verts_size();

    cl_mem src,dst,weights,constraints,indices,tmp_res,bounds,final_res;

    src = m_Memobjs[0];			//Vertices
    constraints = m_Memobjs[1];	//Constraints
    indices = m_Memobjs[2];		//Indices
    weights = m_Memobjs[3];		//Weights
    dst = m_Memobjs[4];			//Output
    tmp_res = m_Memobjs[5];
    bounds = m_Memobjs[6];
    final_res = m_Memobjs[7];
    
    //set kernel arguments
    cl_int err = clSetKernelArg(m_Kernel, 0, sizeof(cl_mem), (void *)&src );
    err |= clSetKernelArg(m_Kernel, 1, sizeof(cl_mem), (void *)&constraints );
    err |= clSetKernelArg(m_Kernel, 2, sizeof(cl_mem), (void *)&indices );
    err |= clSetKernelArg(m_Kernel, 3, sizeof(cl_mem), (void *)&dst );
    err |= clSetKernelArg(m_Kernel, 4, sizeof(cl_int), &n );

    if( m_Options.weights == 1 )	//conformal weights
        err |= clSetKernelArg(m_Kernel, 5, sizeof(cl_mem), (void *)&weights );
    else
        err |= clSetKernelArg(m_Kernel, 5, sizeof(cl_mem), (void *)&bounds );

    err |= clSetKernelArg(m_KernelSPConvergence, 0, sizeof(cl_mem), (void *)&src );
    err |= clSetKernelArg(m_KernelSPConvergence, 1, sizeof(cl_mem), (void *)&constraints );
    err |= clSetKernelArg(m_KernelSPConvergence, 2, sizeof(cl_mem), (void *)&indices );
    err |= clSetKernelArg(m_KernelSPConvergence, 3, sizeof(cl_mem), (void *)&dst );
    err |= clSetKernelArg(m_KernelSPConvergence, 4, sizeof(cl_int), &n );
    err |= clSetKernelArg(m_KernelSPConvergence, 5, sizeof(cl_mem), (void *)&tmp_res );

    if( m_Options.weights == 1 )	//conformal weights
        err |= clSetKernelArg(m_KernelSPConvergence, 6, sizeof(cl_mem), (void *)&weights );
    else
        err |= clSetKernelArg(m_KernelSPConvergence, 6, sizeof(cl_mem), (void *)&bounds );

    err |= clSetKernelArg(m_KernelConvergence, 1, sizeof(cl_mem), (void *)&indices );
    err |= clSetKernelArg(m_KernelConvergence, 2, sizeof(cl_int), &n );
    err |= clSetKernelArg(m_KernelConvergence, 3, sizeof(cl_mem), (void *)&tmp_res );

    if(  m_Options.weights == 1 )
        err |= clSetKernelArg(m_KernelConvergence, 4, sizeof(cl_mem), (void *)&weights );
    else
        err |= clSetKernelArg(m_KernelConvergence, 4, sizeof(cl_mem), (void *)&bounds );

    if( err != CL_SUCCESS )
    {
        std::cerr << "Some error occured during the setting of arguments : " << err << std::endl;
        return false;
    }

    cl_scalar_t last_res = 1000000.0f;
    cl_scalar_t res, linear_res;
    std::size_t one = 1;

    for( std::size_t i = 0; i < m_Options.max_iters; ++i )
    {
        solve_stats.iterations++;

        for( std::size_t j = 1; j < m_Options.max_sp_iters; ++j )
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
                clSetKernelArg(m_Kernel, 3, sizeof(cl_mem), (void *)&dst );
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
                clSetKernelArg(m_Kernel, 3, sizeof(cl_mem), (void *)&dst );

                //Finally compute the residual (Lmax)
                std::size_t s = m_ReductionBlocks * m_ReductionThreads;
    
                //run the max reduction kernels
                while( s > 1 )
                {
                    clSetKernelArg(m_KernelSRes, 0, sizeof(cl_int), &n );
                    clSetKernelArg(m_KernelSRes, 1, sizeof(cl_mem), &tmp_res );
                    clSetKernelArg(m_KernelSRes, 2, sizeof(cl_mem), &final_res);
                    clSetKernelArg(m_KernelSRes, 3, sizeof(cl_scalar_t)*m_ReductionThreads, NULL);

                    std::size_t threads = (s < m_ReductionThreads) ? cl::next_pow2(s) : m_ReductionThreads;
                    std::size_t blocks = (s + threads - 1) / threads;
                 
                    std::size_t rglobal_work_size = threads * blocks;
                    std::size_t rlocal_work_size = threads;

                    err |= clEnqueueNDRangeKernel(m_CmdQueue, m_KernelSRes, 1, NULL, &rglobal_work_size, &rlocal_work_size, 0, NULL, NULL);

                    s = (s + threads - 1) / threads;
                }
              
                if( err != CL_SUCCESS )
                {
                    std::cerr << "Some error occured during the execution of the kernel : " << err << std::endl;
                    return false;
                }

                //Read back the residual value
                err |= clEnqueueReadBuffer(m_CmdQueue, final_res, CL_TRUE, 0, sizeof(cl_scalar_t), &linear_res, 0, NULL, NULL);

                if( linear_res < m_Options.spdelta )
                    break;				
            }
        }

        //Normalize the solution and update the constraints
        clSetKernelArg(m_KernelNormalize, 0, sizeof(cl_mem), (void *)&src);			//vertices
        clSetKernelArg(m_KernelNormalize, 1, sizeof(cl_int), (void *)&n);			//vertices
        clSetKernelArg(m_KernelNormalize, 2, sizeof(cl_mem), (void *)&constraints);	//constraints
        clSetKernelArg(m_KernelNormalize, 3, sizeof(cl_mem), (void *)&dst);			//vertices

        err = clEnqueueNDRangeKernel(m_CmdQueue, m_KernelNormalize, 1, NULL, &m_GlobalSize, &m_LocalSize, 0, NULL, NULL);

        //swap buffers
        std::swap(dst,src);

        //Calculate the residual for each vertex
        clSetKernelArg(m_KernelConvergence, 0, sizeof(cl_mem), (void *)&src );

        err = clEnqueueNDRangeKernel(m_CmdQueue, m_KernelConvergence, 1, NULL, &m_GlobalSize, &m_LocalSize, 0, NULL, NULL);

        //Finally compute the residual (L2)
        clSetKernelArg(m_KernelRes, 0, sizeof(cl_int), &n );
        clSetKernelArg(m_KernelRes, 1, sizeof(cl_mem), &tmp_res );
        clSetKernelArg(m_KernelRes, 2, sizeof(cl_mem), &final_res);

        err = clEnqueueNDRangeKernel(m_CmdQueue, m_KernelRes, 1, NULL, &one, &one, 0, NULL, NULL);

        //Read back the residual value
        err = clEnqueueReadBuffer(m_CmdQueue, final_res, CL_TRUE, 0, sizeof(cl_scalar_t), &res, 0, NULL, NULL);

        std::cout << i << " - res : " << res << '\n';

        if( last_res < res || res < m_Options.target_residual )
            break;

        last_res = res;
    }

    //Read the final coordinates
    err = clEnqueueReadBuffer(m_CmdQueue, src, CL_TRUE, 0, n*sizeof(cl_scalar4_t), m_SVertices, 0, NULL, NULL);

    solve_stats.elapsed_ms = timeGetTime() - solve_stats.elapsed_ms;


    for( std::size_t i = 0; i < n; ++i )
    {
        m_Output.add_vertex( vertex_ptr_t( new vertex_t(
            m_SVertices[i].s[0],m_SVertices[i].s[1],m_SVertices[i].s[2] ) ) );
    }

    for( mesh_t::face_iterator_t f = m_Input.faces_begin(); f != m_Input.faces_end(); ++f )
    {
        std::vector< vertex_ptr_t > verts;
        std::vector< vertex_ptr_t > indices((*f)->verts_begin(),(*f)->verts_end());

        for( std::size_t i = 0 ; i < indices.size(); ++i )
            verts.push_back( m_Output.verts()[ indices[i]->get_id() ] );

        m_Output.add_face( face_ptr_t( new face_t(verts.begin(),verts.end()) ) );
    }

    solve_stats.residual = res;

    return true;
}

void spheremap::SolverCL::compute_weights()
{
    //Compute Conformal or Barycentric weights for all the edges
    foreach( mesh_t::edge_ptr_t e, m_Input.edges() )
    {
        std::vector< vertex_ptr_t > verts(e->verts_begin(),e->verts_end());

        //Find the two oposing vertices
        foreach( mesh_t::face_ptr_t ef, e->faces() )
        {
            foreach( mesh_t::vertex_ptr_t efv, ef->verts() )
            {
                if( efv != verts[0] && efv != verts[1] )
                    verts.push_back(efv);
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
    foreach( mesh_t::vertex_ptr_t v, m_Input.verts() )
    {
        float sum = 0.0;
        for( std::size_t j = 0; j < v->Weights.size(); ++j )
            sum += v->Weights[j];

        for( std::size_t j = 0; j < v->Weights.size(); ++j )
            v->Weights[j] /= sum;
    }
}
