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

#ifndef OPENCL_BLAS_H
#define OPENCL_BLAS_H

#define __CL_ENABLE_EXCEPTIONS
#include <CL/cl.hpp>
#include "config.h"

#if CL_USE_DOUBLE
typedef cl_double cl_scalar_t;
#else
typedef cl_float cl_scalar_t;
#endif

namespace arch
{

namespace blas
{

enum status_t
{
	STATUS_SUCCESS,
	STATUS_ALLOC_FAILED,
	STATUS_EXECUTION_FAILED
};

struct handler 
{
	handler(cl::CommandQueue &queue);

	cl::Kernel dot;
	cl::Kernel axpy;
	cl::Kernel scal;
	cl::Kernel sum;
    cl::Kernel avg;
    cl::Kernel xmy;
	cl::Buffer TempValue;
    cl::Event event;
	cl::CommandQueue &BlasQueue;
    cl::Context BlasContext;
    cl::Device BlasDevice;
    cl::Buffer m_temp_buffer;
    std::size_t m_temp_buffer_size;
	cl_int MaxWorksize;
};

status_t dot(
    handler &h,
	const int n,
	const cl::Buffer &x,
    const int incx,
    cl::Buffer &y,
    int incy,
    cl_scalar_t *result);

status_t axpy(
	handler &h,
    const int n,
    const cl_scalar_t alpha,
    const cl::Buffer &x,
    const int incx,
    cl::Buffer &y,
    const int incy);

status_t xmy(
    handler &h,
    const int n,
    const cl::Buffer &x,
    const cl::Buffer &y,
    const cl::Buffer &z);

status_t scal(
    handler &h,
    const int n,
    const cl_scalar_t alpha,
    cl::Buffer &x,
    const int incx);

status_t copy(
     handler &h,
	 const int n,
	 cl::Buffer &x,
	 const int incx,
	 cl::Buffer &y,
	 const int incy);

status_t sum(
	 handler &h,
	 const int n,
	 cl::Buffer &x,
	 cl::Buffer &tmp,
   cl_scalar_t *result);

status_t avg(
    handler &h,
    const int n,
    cl::Buffer &x,
    cl_scalar_t *result);

status_t amax(
    handler &h,
    const int n,
    cl::Buffer &x,
    cl_scalar_t *result);

status_t amin(
    handler &h,
    const int n,
    cl::Buffer &x,
    cl_scalar_t *result);

void print(
	 handler &h,
	 const int n,
	 cl::Buffer &x,
   bool vertical = false);

void profile(int sel_platform, int sel_device);

double stable_sum( const int n, const double *vec );
double stable_sum( const int n, const float *vec );

}

}


#endif
