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

#ifndef USE_DOUBLE
#define USE_DOUBLE 0
#endif

#if USE_DOUBLE
#ifdef cl_khr_fp64
    #pragma OPENCL EXTENSION cl_khr_fp64 : enable
#elif defined(cl_amd_fp64)
    #pragma OPENCL EXTENSION cl_amd_fp64 : enable
#else 
    #error "Double precision not supported by OpenCL implementation"
#endif

typedef double scalar_t;
typedef double2 scalar2_t;
typedef double4 scalar4_t;
typedef double8 scalar8_t;
#else
typedef float scalar_t;
typedef float2 scalar2_t;
typedef float4 scalar4_t;
typedef float8 scalar8_t;
#endif

#define M_SQRT1_3 0.5773502691896258

kernel void laplace_kernel(
     const int n,
     const global scalar2_t *pnts,
     const global int *ia,
     const global int *ja,
     global scalar2_t *result)
{
  int i = get_global_id(0);
  if( i < n ) 
  {
    scalar2_t c = (scalar2_t)(0.0,0.0);
    for( int j = ia[i]; j < ia[i+1]; ++j ) {
      c += pnts[ ja[ j ] ];
    }
    c /= (scalar_t)(ia[i+1] - ia[i]);
    result[i] = c;
  }
}

kernel void untangle_kernel(
     const int n,
     const global scalar2_t *pnts,
     const global int *ia,
     const global int *ja,
     const double sigma,
     global scalar2_t *result)
{
  int i = get_global_id(0);
  if( i < n ) 
  {
    scalar2_t c = (scalar2_t)(0.0,0.0);
    int start = ia[i];
    int end = ia[i+1];
    scalar_t w = 1.0;
  
    scalar2_t v0 = pnts[ start ];
    scalar2_t v1 = pnts[ end ];
    scalar2_t da;
    da.x = v1.y - v0.y;
    da.y = v0.x - v1.x;
    scalar2_t c = w * v0 + sigma * da;

    for( j = start + 1; j < end; ++j ) {
      v1 = pnts[ ja[j] ];
      scalar2_t da;
      da.x = v0.y - v1.y;
      da.y = v1.x - v0.x;
      c += w * v0 + sigma * da;
      v0 = v1;
    }

    c /= ((1.0 + sigma) * (scalar_t)(ia[i+1] - ia[i]));
    result[i] = c;
  }
}

