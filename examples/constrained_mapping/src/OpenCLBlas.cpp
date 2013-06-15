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

#include "OpenCLBlas.h"
#include "UtilsCL.h"
#include <iostream>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real_distribution.hpp>

using namespace arch::clutils;

#define CL_STRINGIFY(A) #A

const char blas_source[] = 
    CL_STRINGIFY(
#if USE_DOUBLE\n
#ifdef cl_khr_fp64\n
#pragma OPENCL EXTENSION cl_khr_fp64 : enable\n
#elif defined(cl_amd_fp64)\n
#pragma OPENCL EXTENSION cl_amd_fp64 : enable\n
#else \n
#error "Double precision not supported by OpenCL implementation"\n
#endif\n
    typedef double scalar_t;\n
#else\n
    typedef float scalar_t;\n
#endif\n

 kernel void blasdot(
 const int n,
 global const scalar_t *x,
 const int incx,
 global const scalar_t *y,
 const int incy,
 local scalar_t *work,
 global scalar_t *result)
{
    //for the recuction the size of the workgroup must be a power of 2
    unsigned int tid = get_local_id(0);
    unsigned int i = get_group_id(0)*(2*get_local_size(0)) + get_local_id(0);

    work[tid] = (i < n) ? x[i]*y[i] : 0;
    if( (i + get_local_size(0)) < n ) 
        work[tid] += x[i + get_local_size(0)]*y[i + get_local_size(0)];

    barrier(CLK_LOCAL_MEM_FENCE);

    //do reduction in shared mem
    for( unsigned int s = get_local_size(0)/2; s > 0; s >>=1 ) 
    {
        if( tid < s )
            work[tid] += work[tid + s];
    
        barrier(CLK_LOCAL_MEM_FENCE);
    }
   
    if(tid == 0) result[get_group_id(0)] = work[0];
}

kernel void blassum(
    const int n, 
    const __global scalar_t *x, 
    local scalar_t *work,
     global scalar_t *y)
{
    //for the recuction the size of the workgroup must be a power of 2
    unsigned int tid = get_local_id(0);
    unsigned int i = get_group_id(0)*(2*get_local_size(0)) + get_local_id(0);

    work[tid] = (i < n) ? x[i] : 0;
    if( (i + get_local_size(0)) < n ) 
        work[tid] += x[i + get_local_size(0)];

    barrier(CLK_LOCAL_MEM_FENCE);

    //do reduction in shared mem
    for( unsigned int s = get_local_size(0)/2; s > 0; s >>=1 ) 
    {
        if( tid < s )
            work[tid] += work[tid + s];

        barrier(CLK_LOCAL_MEM_FENCE);
    }

    if(tid == 0) y[get_group_id(0)] = work[0];
}

kernel void blasavg(
    const int n, 
    const __global scalar_t *x, 
    local scalar_t *work,
    global scalar_t *y)
{
    //for the recuction the size of the workgroup must be a power of 2
    unsigned int tid = get_local_id(0);
    unsigned int i = get_group_id(0)*(2*get_local_size(0)) + get_local_id(0);

    work[tid] = (i < n) ? x[i] : 0;
    if( (i + get_local_size(0)) < n ) 
        work[tid] += x[i + get_local_size(0)];

    barrier(CLK_LOCAL_MEM_FENCE);

    //do reduction in shared mem
    for( unsigned int s = get_local_size(0)/2; s > 0; s >>=1 ) 
    {
        if( tid < s )
            work[tid] += work[tid + s];

        barrier(CLK_LOCAL_MEM_FENCE);
    }

    if(tid == 0) y[get_group_id(0)] = (work[0] / n);
}

 kernel void blasdot_nv(
 const int n,
 global const scalar_t *x,
 const int incx,
 global const scalar_t *y,
 const int incy,
 local volatile scalar_t *work,
 global scalar_t *result)
{
    //for the recuction the size of the workgroup must be a power of 2
    unsigned int tid = get_local_id(0);
    unsigned int i = get_group_id(0)*(2*BLAS_BLOCK_SIZE) + get_local_id(0);
    double priv_acc;

    priv_acc = (i < n) ? x[i]*y[i] : 0;
    i += BLAS_BLOCK_SIZE;
    if( i < n ) 
        priv_acc += x[i]*y[i];

    work[tid] = priv_acc;
    barrier(CLK_LOCAL_MEM_FENCE);
        
\n#if BLAS_BLOCK_SIZE >= 512                      \n   
    if( tid < 256 ) { priv_acc += work[tid + 256]; work[tid] = priv_acc;}\n
    barrier(CLK_LOCAL_MEM_FENCE);                 \n
#endif                                            \n

\n#if BLAS_BLOCK_SIZE >= 256                      \n   
    if( tid < 128 ) { priv_acc += work[tid + 128]; work[tid] = priv_acc;}\n
    barrier(CLK_LOCAL_MEM_FENCE);                 \n
#endif                                            \n
    
\n#if BLAS_BLOCK_SIZE >= 128                      \n
    if( tid < 64 ) { priv_acc += work[tid + 64]; work[tid] = priv_acc;}\n
    barrier(CLK_LOCAL_MEM_FENCE);                 \n
#endif                                            \n

    if( tid < 32 ) {
        priv_acc += work[tid + 32]; work[tid] = priv_acc;
        priv_acc += work[tid + 16]; work[tid] = priv_acc;
        priv_acc += work[tid +  8]; work[tid] = priv_acc;
        priv_acc += work[tid +  4]; work[tid] = priv_acc;
        priv_acc += work[tid +  2]; work[tid] = priv_acc;
        priv_acc += work[tid +  1]; work[tid] = priv_acc;
    }

    if(tid == 0) result[get_group_id(0)] = priv_acc;
}

 kernel void blassum_nv(
    const int n, 
    const global scalar_t *x, 
    local volatile scalar_t *work,
    global scalar_t *result)
{
    //for the recuction the size of the workgroup must be a power of 2
    unsigned int tid = get_local_id(0);
    unsigned int i = get_group_id(0)*(2*BLAS_BLOCK_SIZE) + get_local_id(0);
    double priv_acc;

    priv_acc = (i < n) ? x[i] : 0;
    if( (i + BLAS_BLOCK_SIZE) < n ) 
        priv_acc += x[i + BLAS_BLOCK_SIZE];

    work[tid] = priv_acc;
    barrier(CLK_LOCAL_MEM_FENCE);
    
\n#if BLAS_BLOCK_SIZE >= 512                      \n   
    if( tid < 256 ) { priv_acc += work[tid + 256]; work[tid] = priv_acc;}\n
    barrier(CLK_LOCAL_MEM_FENCE);                 \n
#endif                                            \n

\n#if BLAS_BLOCK_SIZE >= 256                      \n   
    if( tid < 128 ) { priv_acc += work[tid + 128]; work[tid] = priv_acc;}\n
    barrier(CLK_LOCAL_MEM_FENCE);                 \n
#endif                                            \n
    
\n#if BLAS_BLOCK_SIZE >= 128                      \n
    if( tid < 64 ) { priv_acc += work[tid + 64]; work[tid] = priv_acc;}\n
    barrier(CLK_LOCAL_MEM_FENCE);                 \n
#endif                                            \n

    if( tid < 32 ) {
        priv_acc += work[tid + 32]; work[tid] = priv_acc;
        priv_acc += work[tid + 16]; work[tid] = priv_acc;
        priv_acc += work[tid +  8]; work[tid] = priv_acc;
        priv_acc += work[tid +  4]; work[tid] = priv_acc;
        priv_acc += work[tid +  2]; work[tid] = priv_acc;
        priv_acc += work[tid +  1]; work[tid] = priv_acc;
    }

    if(tid == 0) result[get_group_id(0)] = priv_acc;
}

 kernel void blasavg_nv(
    const int n, 
    const global scalar_t *x, 
    local volatile scalar_t *work,
    global scalar_t *y)
{
    //for the recuction the size of the workgroup must be a power of 2
    unsigned int tid = get_local_id(0);
    unsigned int i = get_group_id(0)*(2*BLAS_BLOCK_SIZE) + get_local_id(0);
    double priv_acc;

    priv_acc = (i < n) ? x[i] : 0;
    if( (i + BLAS_BLOCK_SIZE) < n ) 
        priv_acc += x[i + BLAS_BLOCK_SIZE];
       
    work[tid] = priv_acc;
    barrier(CLK_LOCAL_MEM_FENCE);
    
\n#if BLAS_BLOCK_SIZE >= 512                      \n   
    if( tid < 256 ) { priv_acc += work[tid + 256]; work[tid] = priv_acc;}\n
    barrier(CLK_LOCAL_MEM_FENCE);                 \n
#endif                                            \n

\n#if BLAS_BLOCK_SIZE >= 256                      \n   
    if( tid < 128 ) { priv_acc += work[tid + 128]; work[tid] = priv_acc;}\n
    barrier(CLK_LOCAL_MEM_FENCE);                 \n
#endif                                            \n
    
\n#if BLAS_BLOCK_SIZE >= 128                      \n
    if( tid < 64 ) { priv_acc += work[tid + 64]; work[tid] = priv_acc;}\n
    barrier(CLK_LOCAL_MEM_FENCE);                 \n
#endif                                            \n

    if( tid < 32 ) {
        priv_acc += work[tid + 32]; work[tid] = priv_acc;
        priv_acc += work[tid + 16]; work[tid] = priv_acc;
        priv_acc += work[tid +  8]; work[tid] = priv_acc;
        priv_acc += work[tid +  4]; work[tid] = priv_acc;
        priv_acc += work[tid +  2]; work[tid] = priv_acc;
        priv_acc += work[tid +  1]; work[tid] = priv_acc;
    }

    if(tid == 0) y[get_group_id(0)] = (priv_acc / n);
}

kernel void blasaxpy_ni(
    const int n,
    const scalar_t alpha,
    global const scalar_t *x,
    global scalar_t *y)
{
    int tid = get_global_id(0);
    if( tid < n ) y[tid] = mad( alpha, x[tid], y[tid] );
}

kernel void blasscal_ni(
    const int n,
    const scalar_t alpha,
    global scalar_t *x)
{
    int tid = get_global_id(0);
    if( tid < n ) x[tid] = x[tid] * alpha;
}

kernel void blasxmy(
    const int n,
    global const scalar_t *x,
    global const scalar_t *y,
    global scalar_t *z)
{
    int tid = get_global_id(0);
    if( tid < n ) z[tid] = x[tid] - y[tid];
}

kernel void blasaxpy(
    const int n,
    const scalar_t alpha,
    global const scalar_t *x,
    const int incx,
    global scalar_t *y,
    const int incy)
{
    int tid = get_global_id(0);
    if( tid < n ) y[tid*incy] = mad( alpha, x[tid*incx], y[tid*incy] );
}

kernel void blasscal(
    const int n,
    const scalar_t alpha,
    global scalar_t *x,
    const int incx)
{
    int tid = get_global_id(0);
    if( tid < n ) x[tid*incx] *= alpha;
}
);

#if 0
#include <immintrin.h>

double stable_sum( const int n, const double *vec )
{
     double result = 0.0;
    __m256d sum = _mm256_set_pd(0.0,0.0,0.0,0.0);
    __m256d c   = _mm256_set_pd(0.0,0.0,0.0,0.0);

    const int round4 = (n / 4) * 4;
    for( int i = 0; i < round4; i += 4 ) {
        __m256d y;
        __m256d t;

        __m256d v = _mm256_load_pd( &vec[i] );
        y = _mm256_sub_pd( v, c );
        t = _mm256_add_pd( sum, y );
        c = _mm256_sub_pd( t, sum);
        c = _mm256_sub_pd( c, y );
        sum = t;
    }

    result = sum.m256d_f64[0] + sum.m256d_f64[1] + sum.m256d_f64[2] + sum.m256d_f64[3];
            
    for( int i = round4; i < n; ++i ) 
        result += vec[i];

    return result;
}

#else
double arch::blas::stable_sum( const int n, const double *vec )
{
    double result = 0.0;
    double sum[4] = {0.0};
    double c[4]   = {0.0};
    const int round4 = (n / 4) * 4;
    for( int i = 0; i < round4; i += 4 ) {
        double y[4];
        double t[4]; 

        y[0] = vec[i + 0] - c[0];
        y[1] = vec[i + 1] - c[1];
        y[2] = vec[i + 2] - c[2];
        y[3] = vec[i + 3] - c[3];

        t[0] = sum[0] + y[0];
        t[1] = sum[1] + y[1];
        t[2] = sum[2] + y[2];
        t[3] = sum[3] + y[3];

        c[0] = (t[0] - sum[0]) - y[0];
        c[1] = (t[1] - sum[1]) - y[1];
        c[2] = (t[2] - sum[2]) - y[2];
        c[3] = (t[3] - sum[3]) - y[3];

        sum[0] = t[0];
        sum[1] = t[1];
        sum[2] = t[2];
        sum[3] = t[3];
    }
    result = sum[0] + sum[1] + sum[2] + sum[3];
            
    for( int i = round4; i < n; ++i ) 
        result += vec[i];

    return result;
}

double arch::blas::stable_sum( const int n, const float *vec )
{
    float result = 0.0;
    float sum[4] = {0.0};
    float c[4]   = {0.0};
    const int round4 = (n / 4) * 4;
    for( int i = 0; i < round4; i += 4 ) {
        float y[4];
        float t[4]; 

        y[0] = vec[i + 0] - c[0];
        y[1] = vec[i + 1] - c[1];
        y[2] = vec[i + 2] - c[2];
        y[3] = vec[i + 3] - c[3];

        t[0] = sum[0] + y[0];
        t[1] = sum[1] + y[1];
        t[2] = sum[2] + y[2];
        t[3] = sum[3] + y[3];

        c[0] = (t[0] - sum[0]) - y[0];
        c[1] = (t[1] - sum[1]) - y[1];
        c[2] = (t[2] - sum[2]) - y[2];
        c[3] = (t[3] - sum[3]) - y[3];

        sum[0] = t[0];
        sum[1] = t[1];
        sum[2] = t[2];
        sum[3] = t[3];
    }
    result = sum[0] + sum[1] + sum[2] + sum[3];
            
    for( int i = round4; i < n; ++i ) 
        result += vec[i];

    return result;
}
#endif

arch::blas::status_t arch::blas::avg(
    arch::blas::handler &h,
    const int n,
    cl::Buffer &x,
    cl_scalar_t *result)
{
    try
    {
        cl_int reduction_threads = h.MaxWorksize;
        cl_int work_size = 2*reduction_threads;
        cl_int threads = reduction_threads;
        cl_int blocks = (n + (work_size - 1)) / (work_size);
      
        if( h.m_temp_buffer_size < blocks ) {
            h.m_temp_buffer = cl::Buffer(h.BlasContext, CL_MEM_READ_WRITE, sizeof(double)*blocks );
            h.m_temp_buffer_size = blocks;
        }
        
        h.avg.setArg(0,n);
        h.avg.setArg(1,x);
        h.avg.setArg(2,sizeof(cl_scalar_t)*reduction_threads, NULL );
        h.avg.setArg(3,h.m_temp_buffer);

        h.BlasQueue.enqueueNDRangeKernel(
            h.avg,
            cl::NullRange,
            cl::NDRange(threads*blocks),
            cl::NDRange(threads),
            NULL,
            &h.event
            );

        h.event.wait();

        if( blocks < 4096 ) {
            std::vector< cl_double > temp( blocks );
            h.BlasQueue.enqueueReadBuffer( h.m_temp_buffer, GL_TRUE, 0, sizeof(cl_double) * blocks, &temp[0] );
            *result = stable_sum( blocks, &temp[0] );
        } else {
            sum(h, blocks, h.m_temp_buffer, h.m_temp_buffer, result);
        }

        //sum(h, blocks, d_temp, d_temp, result);
    } catch ( cl::Error &err ) {
        std::cerr << "arch::blas::avg :: Error in : " << err.what() << " code : " << err.err() << "\n";
        return arch::blas::STATUS_EXECUTION_FAILED;
    }

    return arch::blas::STATUS_SUCCESS;
}

arch::blas::status_t arch::blas::dot(
    arch::blas::handler &h,
    const int n,
    const cl::Buffer &x,
    const int incx,
    cl::Buffer &y,
    const int incy,
    cl_scalar_t *result)
{
    try
    {
        cl_int reduction_threads = h.MaxWorksize;
        cl_int work_size = 2*reduction_threads;
        cl_int threads = reduction_threads;
        cl_int blocks = (n + (work_size - 1)) / (work_size);
      
        if( h.m_temp_buffer_size < blocks ) {
            h.m_temp_buffer = cl::Buffer(h.BlasContext, CL_MEM_READ_WRITE, sizeof(double)*blocks );
            h.m_temp_buffer_size = blocks;
        }

        h.dot.setArg(0,n);
        h.dot.setArg(1,x);
        h.dot.setArg(2,incx);
        h.dot.setArg(3,y);
        h.dot.setArg(4,incy);
        h.dot.setArg(5,sizeof(cl_scalar_t)*reduction_threads, NULL );
        h.dot.setArg(6,h.m_temp_buffer);

        h.BlasQueue.enqueueNDRangeKernel(
            h.dot,
            cl::NullRange,
            cl::NDRange(threads*blocks),
            cl::NDRange(threads),
            NULL,
            &h.event
            );

        //h.event.wait();

        //Do Kahan on the CPU since it is faster for small data sets
        if( blocks < 4096 ) {
            std::vector< cl_double > temp( blocks );
            h.BlasQueue.enqueueReadBuffer( h.m_temp_buffer, GL_TRUE, 0, sizeof(cl_double) * blocks, &temp[0] );
            *result = stable_sum( blocks, &temp[0] );
        } else {
            sum(h, blocks, h.m_temp_buffer, h.m_temp_buffer, result);
        }
    } catch ( cl::Error &err ) {
        std::cerr << "arch::blas::sdot :: Error in : " << err.what() << " code : " << err.err() << "\n";
        return arch::blas::STATUS_EXECUTION_FAILED;
    }

    return arch::blas::STATUS_SUCCESS;
}

arch::blas::status_t arch::blas::axpy(
    arch::blas::handler &h,
    const int n,
    const cl_scalar_t alpha,
    const cl::Buffer &x,
    const int incx,
    cl::Buffer &y,
    const int incy)
{

    if( alpha == 0.0 ) return arch::blas::STATUS_SUCCESS;

    try
    {
        h.axpy.setArg(0,n);
        h.axpy.setArg(1,alpha);
        h.axpy.setArg(2,x);
        h.axpy.setArg(3,y);
     
        int nsize = round_up(h.MaxWorksize,n);

        h.BlasQueue.enqueueNDRangeKernel(
            h.axpy,
            cl::NullRange,
            cl::NDRange(nsize),
            cl::NDRange(h.MaxWorksize),
            NULL,
            &h.event
            );

        h.event.wait();

    } catch ( cl::Error &err ) {
        std::cerr << "arch::blas::saxpy :: Error in : " << err.what() << " code : " << err.err() << " error = " << arch::clutils::error(err.err() ) << "\n";
        return arch::blas::STATUS_EXECUTION_FAILED;
    }

    return arch::blas::STATUS_SUCCESS;
}

arch::blas::status_t arch::blas::xmy(
    handler &h,
    const int n,
    const cl::Buffer &x,
    const cl::Buffer &y,
    const cl::Buffer &z)
{
    try
    {
        h.xmy.setArg(0,n);
        h.xmy.setArg(1,x);
        h.xmy.setArg(2,y);
        h.xmy.setArg(3,z);
      
        int nsize = round_up(h.MaxWorksize,n);

        h.BlasQueue.enqueueNDRangeKernel(
            h.xmy,
            cl::NullRange,
            cl::NDRange(nsize),
            cl::NDRange(h.MaxWorksize),
            NULL,
            &h.event
            );

        h.event.wait();

    } catch ( cl::Error &err ) {
        std::cerr << "arch::blas::xmy :: Error in : " << err.what() << " code : " << err.err() << " error = " << arch::clutils::error(err.err() ) << "\n";
        return arch::blas::STATUS_EXECUTION_FAILED;
    }

    return arch::blas::STATUS_SUCCESS;
}

arch::blas::status_t arch::blas::scal(
    arch::blas::handler &h,
    const int n,
    const cl_scalar_t alpha,
    cl::Buffer &x,
    const int incx)
{
    if( alpha == 1.0 ) return arch::blas::STATUS_SUCCESS;

    try
    {
        h.scal.setArg(0,n);
        h.scal.setArg(1,alpha);
        h.scal.setArg(2,x);
      
        int nsize = round_up(h.MaxWorksize,n);

        h.BlasQueue.enqueueNDRangeKernel(
            h.scal,
            cl::NullRange,
            cl::NDRange(nsize),
            cl::NDRange(h.MaxWorksize),
            NULL,
            &h.event
            );

        h.event.wait();

    } catch ( cl::Error &err ) {
        std::cerr << "arch::blas::sscal :: Error in : " << err.what() << " code : " << err.err() << " error = " << arch::clutils::error(err.err() ) << "\n";
        return arch::blas::STATUS_EXECUTION_FAILED;
    }

    return arch::blas::STATUS_SUCCESS;
}

arch::blas::status_t arch::blas::copy( 
    arch::blas::handler &h,
    const int n,
    cl::Buffer &x,
    const int incx,
    cl::Buffer &y,
    const int incy)
{
    try
    {
        h.BlasQueue.enqueueCopyBuffer(x,y,0,0,n*sizeof(cl_scalar_t) );
    } catch ( cl::Error &err ) {
        std::cerr << "arch::blas::scopy :: Error in : " << err.what() << " code : " << err.err() << " error = " << arch::clutils::error(err.err() ) << "\n";
        return arch::blas::STATUS_EXECUTION_FAILED;
    }

    return arch::blas::STATUS_SUCCESS;
}

arch::blas::status_t arch::blas::sum( 
    arch::blas::handler &h,
    const int n,
    cl::Buffer &x,
    cl::Buffer &tmp,
    cl_scalar_t *result)
{
    try
    {
        cl_int reduction_threads = h.MaxWorksize;
        cl_int work_size = 2*reduction_threads;
        cl_int s = n;

        //run the max reduction kernels
        while( s > 1 )
        {
            h.sum.setArg(0,s);
            h.sum.setArg(1,x);
           
            cl_int threads = reduction_threads;
            cl_int blocks = (s + (work_size - 1)) / (work_size);
            h.sum.setArg(2, sizeof(cl_scalar_t)*threads, NULL );
            h.sum.setArg(3,tmp );

            h.BlasQueue.enqueueNDRangeKernel(
                h.sum,
                cl::NullRange,
                cl::NDRange(threads * blocks),
                cl::NDRange(threads),
                NULL,
                &h.event
                );

            h.event.wait();

            s = blocks;

            std::swap(x,tmp);
        }

        h.BlasQueue.enqueueReadBuffer(
            x, CL_TRUE, 0, sizeof(cl_scalar_t), result );

    } catch ( cl::Error &err ) {
        std::cerr << "arch::blas::sum :: Error in : " << err.what() << " code : " << err.err() << " error = " << arch::clutils::error(err.err() ) << "\n";
        return arch::blas::STATUS_EXECUTION_FAILED;
    }

    return arch::blas::STATUS_SUCCESS;
}

arch::blas::status_t arch::blas::amax(
    handler &h,
    const int n,
    cl::Buffer &x,
    cl_scalar_t *result)
{
    std::vector<cl_scalar_t> h_x(n);
    h.BlasQueue.enqueueReadBuffer(x, CL_TRUE, 0, h_x.size()*sizeof(cl_scalar_t), &h_x[0] );
    *result = fabs(h_x[0]);
    for( std::size_t i = 0; i < h_x.size(); ++i )
        *result = max( fabs(h_x[i]), *result );

    return arch::blas::STATUS_SUCCESS;
}

arch::blas::status_t arch::blas::amin(
    handler &h,
    const int n,
    cl::Buffer &x,
    cl_scalar_t *result)
{
    std::vector<cl_scalar_t> h_x(n);
    h.BlasQueue.enqueueReadBuffer(x, CL_TRUE, 0, h_x.size()*sizeof(cl_scalar_t), &h_x[0] );
    *result = fabs(h_x[0]);
    for( std::size_t i = 0; i < h_x.size(); ++i ) {
        *result = min( fabs(h_x[i]), *result );
    }

    return arch::blas::STATUS_SUCCESS;
}

void arch::blas::print(
    arch::blas::handler &h,
    const int n,
    cl::Buffer &x,
    bool vertical)
{
    std::vector<cl_scalar_t> debug(n);
    h.BlasQueue.enqueueReadBuffer(x, CL_TRUE, 0, debug.size()*sizeof(cl_scalar_t), &debug[0] );

    for( std::size_t i = 0; i < debug.size(); ++i ) {
        std::cout << debug[i] << ((vertical || i == n-1) ? "\n" : ",");
    }
}

arch::blas::handler::handler(cl::CommandQueue &queue) : BlasQueue(queue)
{
    using namespace cl;
    using std::min;

    try
    {
        BlasContext = BlasQueue.getInfo<CL_QUEUE_CONTEXT>();
        BlasDevice = BlasQueue.getInfo<CL_QUEUE_DEVICE>();		

        Program::Sources source(1, std::make_pair(blas_source, sizeof(blas_source)+1));
        Program program = Program(BlasContext, source);

        std::vector< Device > devices;
        devices.push_back(BlasDevice);
        //MaxWorksize = 256;
        MaxWorksize = min(256, int(BlasDevice.getInfo<CL_DEVICE_MAX_WORK_GROUP_SIZE>()));

        char options[255];
        try{
            sprintf(options, "-cl-mad-enable -cl-no-signed-zeros -DUSE_DOUBLE=%d -DBLAS_BLOCK_SIZE=%d",CL_USE_DOUBLE,MaxWorksize);
            program.build(devices,options);
        } catch ( Error &err ) {
            std::cerr << program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(BlasDevice) << "\n";
        }

        if( BlasDevice.getInfo<CL_DEVICE_TYPE>() == CL_DEVICE_TYPE_GPU ) {
            dot = Kernel(program, "blasdot_nv");
            sum = Kernel(program, "blassum_nv");
            avg = Kernel(program, "blasavg_nv");
        } else {
            dot = Kernel(program, "blasdot");
            sum = Kernel(program, "blassum");
            avg = Kernel(program, "blasavg");
        }

        axpy = Kernel(program, "blasaxpy_ni");
        scal = Kernel(program, "blasscal_ni");
        xmy  = Kernel(program, "blasxmy");

        TempValue = cl::Buffer(BlasContext, CL_MEM_WRITE_ONLY, sizeof(cl_scalar_t) );
        m_temp_buffer_size = 0;
    } catch ( Error &err )
    {
        std::cerr << "blas_handler() " << "OpenCL error in " << err.what() << ", code = " << err.err() << "\n";
    }
}

bool equal_vecs(arch::blas::handler &hd, cl::Buffer &x, const std::vector<double> &y, double rel_tol = 1e-6)
{
    using std::min;
    using std::max;
    using std::abs;

    std::vector<double> hx(y.size());
    hd.BlasQueue.enqueueReadBuffer(
            x, CL_TRUE, 0, sizeof(double)*y.size(), &hx[0] );

    for( std::size_t i = 0; i < hx.size(); ++i ) 
    {
        if( (abs( hx[i] - y[i] ) / max( abs(hx[i]), abs(y[i]) )) > rel_tol ) 
            return false;
    }
    return true;
}

void test_blas(arch::blas::handler &hd, int n)
{
    using std::min;
    using std::max;
    using std::abs;

    std::cout << "Blas profile and test (n = " << n << ")\n";

    std::vector<double> a(n);
    std::vector<double> b(n);
    std::vector<double> c(n);

    std::cout.precision(12);

    boost::random::mt19937 gen;
    boost::random::uniform_real_distribution<> dist(1.0, 30.0);

    for( int i = 0; i < n; ++i ) 
    {
        a[i] = dist(gen);
        b[i] = dist(gen);
        //std::cout << a[i] << "," << b[i] << "\n";
    }

    cl::Buffer d_a(hd.BlasContext, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(double)*n, &a[0] );
    cl::Buffer d_b(hd.BlasContext, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(double)*n, &b[0] );
    
    arch::blas::axpy(hd, n, 2.0, d_a, 1, d_b, 1);
    for( int i = 0; i < n; ++i ) 
        b[i] += 2.0*a[i];

    std::cout << "blas axpy test [" << (equal_vecs(hd, d_b, b) ? "PASSED" : "FAILED") << "]\n";

    for( int i = 0; i < n; ++i )
        a[i] *= 3.0;
    arch::blas::scal(hd, n, 3.0, d_a, 1);
    std::cout << "blas scal test [" << (equal_vecs(hd, d_a, a) ? "PASSED" : "FAILED") << "]\n";

    double dot_aa = 0.0;
    for( int i = 0; i < n; ++i )
        dot_aa += a[i] * a[i];

    double d_dot_aa;
    arch::blas::dot(hd, n, d_a, 1, d_a, 1, &d_dot_aa);
    double rel_tol = (abs( d_dot_aa - dot_aa ) / max( abs(d_dot_aa), abs(dot_aa) ));
    std::cout << "blas  dot test [" << (rel_tol < 1e-06 ? "PASSED" : "FAILED") << "](" << d_dot_aa << "," << dot_aa << ")(" << rel_tol << ")\n";

    double sum = 0;
    double d_sum = 0;
    for( int i = 0; i < n; ++i )
        sum += b[i];
    arch::blas::sum(hd, n, d_b, d_b, &d_sum);
    rel_tol = (abs( d_sum - sum ) / max( abs(d_sum), abs(d_sum) ));
    std::cout << "blas  sum test [" << (rel_tol < 1e-06 ? "PASSED" : "FAILED") << "](" << d_sum << "," << sum << ")(" << rel_tol << ")\n";

    sum = 0;
    for( int i = 0; i < n; ++i )
        sum += a[i];
    sum /= n;
    arch::blas::avg(hd, n, d_a, &d_sum);
    rel_tol = (abs( d_sum - sum ) / max( abs(d_sum), abs(d_sum) ));
    std::cout << "blas  avg test [" << (rel_tol < 1e-06 ? "PASSED" : "FAILED") << "](" << d_sum << "," << sum << ")(" << rel_tol << ")\n";
}

void arch::blas::profile(int sel_platform, int sel_device)
{
    using namespace cl;

    std::vector<Platform> platforms;
    std::vector<cl::Device> devices;
    Platform::get(&platforms);

    cl_context_properties properties[] = 
    { CL_CONTEXT_PLATFORM, (cl_context_properties)(platforms[sel_platform])(), 0};
    
    Context context(CL_DEVICE_TYPE_ALL, properties); 
    devices = context.getInfo<CL_CONTEXT_DEVICES>();
    
    //create the command queue
    CommandQueue queue(context, devices[sel_device],CL_QUEUE_PROFILING_ENABLE);
    
    handler hd(queue);
    test_blas(hd, 53);
    test_blas(hd, 367);
    test_blas(hd, 997);
    test_blas(hd, 9973);
    test_blas(hd, 100003);
    test_blas(hd, 999983);
}
