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

#ifndef PTHETA
#define PTHETA 1.0   // real value...
#define THETA_IS_ONE 1
#endif

#ifndef PBETA
#define PBETA 0.0
#endif

#define UN_LIMIT 0.0
#define DIR_WEIGHT 1.0
#define UN_WEIGHT 1.0
#define UNTANGLE_FUNC1 0

scalar_t area_theta(scalar_t area)
{
#if THETA_IS_ONE
  return area;
#else
  return powr( area, PTHETA );
#endif
}

scalar_t area_theta_1(scalar_t area)
{
#if THETA_IS_ONE
  return (scalar_t)1.0;
#else
  return powr( area, PTHETA - 1.0 );
#endif
}



#define M_SQRT1_3 0.5773502691896258

kernel void fill_buffer(
       const int n,
       global scalar_t *x,
       const scalar_t pattern)
{
    int tid = get_local_id(0);
    if( tid < n ) 
        x[tid] = pattern;
}

scalar2_t imr_triangle_grad(
       const scalar2_t p0,
       const scalar2_t p1, 
       const scalar2_t p2,
       const scalar_t d42)
{
    scalar4_t a;
    scalar_t h, f, v_sqrt;

    a.even = p1 - p0;
    a.odd  = (2.0*p2 - p1 - p0) * M_SQRT1_3;
    h = (a.s0 * a.s3) - (a.s1 * a.s2);      //determinant
    v_sqrt = sqrt(h*h + d42);
    h += v_sqrt;
    f = dot(a,a);       //frobenius
    h = f / h;
    return (scalar2_t)(
                       2.0*h*(((a.s3             - (a.s2 * M_SQRT1_3))/(2.0*v_sqrt)) - ((a.s0 + (a.s1*M_SQRT1_3)) / f)),
                       2.0*h*((((a.s0*M_SQRT1_3) - a.s1              )/(2.0*v_sqrt)) - ((a.s2 + (a.s3*M_SQRT1_3)) / f)));
}

scalar2_t mips_triangle_grad(
       const scalar2_t p0,
       const scalar2_t p1, 
       const scalar2_t p2,
       const scalar_t cota,
       const scalar_t cotb,
       const scalar_t d42)
{
    scalar4_t a;
    scalar_t h, f, v_sqrt;
    //a.even = p1 - p0 - cota*(p2-p0);
    //a.odd = (p2 - p0)*(cota + cotb);
    a.even = p1 - p0;
    a.odd = (p2 - p1)*cota + (p2 - p0)*cotb;
    h = (a.s0 * a.s3) - (a.s1 * a.s2);      //determinant
    h = max(1e-06,h);  
    v_sqrt = sqrt(h*h + d42);
    //v_sqrt = h;
    h += v_sqrt;
    f = dot(a,a);       //frobenius
    h = f / h;
    
    return (scalar2_t)(
                       2.0*h*((((p1.y-p2.y)*cota + (p1.y-p2.y)*cotb)/(-2.0*v_sqrt)) + ((((p0.x-p1.x) - ((p2.x-p1.x)*cota + (p2.x-p0.x)*cotb)*cotb)) / f)),
                       2.0*h*((((p2.x-p1.x)*cota + (p2.x-p1.x)*cotb)/(-2.0*v_sqrt)) + ((((p0.y-p1.y) - ((p2.y-p1.y)*cota + (p2.y-p0.y)*cotb)*cotb)) / f)));
}

void imr_right_triangle_grad(
       const scalar2_t p0,
       const scalar2_t p1, 
       const scalar2_t p2,
       scalar2_t *grad0,
       scalar2_t *grad1,
       scalar2_t *grad2,
       const scalar_t d42)
{
    //a00, a10, a01, a11
    scalar4_t a;    

    a.lo = p1 - p0;
    a.hi = p2 - p0;
   
    scalar_t frobenius = dot(a,a);
    scalar_t h = (a.s0 * a.s3) - (a.s2 * a.s1);
    scalar_t d_sqrt = sqrt(h*h + d42);
    h = 2*(frobenius / (h + d_sqrt));
    h = h * h;

    //Riva page 66
    grad0[0].s0 += h * ((-a.s0 - a.s2) / frobenius - 0.5*((a.s1-a.s3) / d_sqrt));
    grad0[0].s1 += h * ((-a.s1 - a.s3) / frobenius - 0.5*((a.s2-a.s0) / d_sqrt));
    
    grad1[0].s0 += h * ((a.s0        ) / frobenius - 0.5*((     a.s3) / d_sqrt));
    grad1[0].s1 += h * ((a.s1        ) / frobenius - 0.5*((-a.s2    ) / d_sqrt));

    grad2[0].s0 += h * ((        a.s2) / frobenius - 0.5*((-a.s1    ) / d_sqrt));
    grad2[0].s1 += h * ((        a.s3) / frobenius - 0.5*((     a.s0) / d_sqrt));
    
}

kernel void mips_tris_grad0(
    const int num_of_triangles,
    const scalar_t d42,
    const global scalar2_t *pnts,
    const global int4 *triangles,
    const global scalar4_t *angles,
    global scalar2_t *df_aux)
{
    int tid = get_global_id(0);
    int id = tid*3;
    scalar2_t p0,p1,p2;

    if( tid < num_of_triangles ) 
    {
        int4 f = triangles[tid];
        scalar4_t a = angles[tid]; 
        
        p0 = pnts[ f.s0 ]; 
        p1 = pnts[ f.s1 ]; 
        p2 = pnts[ f.s2 ]; 

        df_aux[ id     ] = mips_triangle_grad(p0,p1,p2,a.s0,a.s1,d42);    //cota,cotb
        df_aux[ id + 1 ] = mips_triangle_grad(p1,p2,p0,a.s1,a.s2,d42);    //cotb,cotc
        df_aux[ id + 2 ] = mips_triangle_grad(p2,p0,p1,a.s2,a.s0,d42);    //cotc,cota
    }
}

kernel void combined_tris_grad0(
    const int num_of_triangles,
    const scalar_t d42,
    const global scalar2_t *pnts,
    const global int4 *triangles,
    const global scalar4_t *angles,
    global scalar2_t *df_aux)
{
    int tid = get_global_id(0);
    int id = tid*3;
    scalar2_t p0,p1,p2;

    if( tid < num_of_triangles ) 
    {
        int4 f = triangles[tid];
        scalar4_t a = angles[tid]; // a.s3 is the area of original triangle
        scalar2_t grad_0, grad_1, grad_2;
        
        p0 = pnts[ f.s0 ]; 
        p1 = pnts[ f.s1 ]; 
        p2 = pnts[ f.s2 ]; 
		
        scalar2_t e0 = p1 - p0;
        scalar2_t e1 = p2 - p1;
        scalar2_t e2 = p2 - p0;
				
        scalar_t det_w = (e0.s0*e2.s1 - e0.s1*e2.s0);
        det_w = max(det_w, 1e-14);
        scalar_t area_energy = (a.s3/det_w) + (det_w/a.s3);
        scalar_t angle_energy = (a.s2*dot(e0,e0) + a.s1*dot(e2,e2) + a.s0*dot(e1,e1))/(2.*det_w); 
        
        grad_0.s0 = -(a.s2*e0.s0  + a.s1*e2.s0 - e1.s1*angle_energy);
        grad_0.s1 = -(a.s2*e0.s1  + a.s1*e2.s1 + e1.s0*angle_energy);
        grad_1.s0 =  (a.s2*e0.s0  - a.s0*e1.s0 - e2.s1*angle_energy);
        grad_1.s1 =  (a.s2*e0.s1  - a.s0*e1.s1 + e2.s0*angle_energy);
        grad_2.s0 =  (a.s0*e1.s0  + a.s1*e2.s0 + e0.s1*angle_energy);
        grad_2.s1 =  (a.s0*e1.s1  + a.s1*e2.s1 - e0.s0*angle_energy);

        area_energy  = (2.0 * area_energy) / det_w;
        angle_energy *= 2.0;
        
        grad_0 = (grad_0 * area_energy);
        grad_1 = (grad_1 * area_energy);
        grad_2 = (grad_2 * area_energy);

        scalar_t factor = (det_w*det_w - a.s3*a.s3) / ((det_w*a.s3)*det_w);

        e0 = (scalar2_t)(-e0.s1,e0.s0);
        e1 = (scalar2_t)(-e1.s1,e1.s0);
        e2 = (scalar2_t)(e2.s1,-e2.s0);

	      df_aux[ id     ] = grad_0 + angle_energy*factor*e1;
        df_aux[ id + 1 ] = grad_1 + angle_energy*factor*e2;
        df_aux[ id + 2 ] = grad_2 + angle_energy*factor*e0;
    }
}

kernel void tris_grad0(
    const int num_of_triangles,
    const scalar_t d42,
    const global scalar2_t *pnts,
    const global int4 *triangles,
    global scalar2_t *df_aux)
{
    int tid = get_global_id(0);
    int id = tid*3;
    scalar2_t p0,p1,p2;

    if( tid < num_of_triangles ) 
    {
        int4 f = triangles[tid];
        p0 = pnts[ f.s0 ]; 
        p1 = pnts[ f.s1 ]; 
        p2 = pnts[ f.s2 ]; 

        df_aux[ id     ] = imr_triangle_grad(p0,p1,p2,d42);
        df_aux[ id + 1 ] = imr_triangle_grad(p1,p2,p0,d42);
        df_aux[ id + 2 ] = imr_triangle_grad(p2,p0,p1,d42);
    }
}

kernel void quads_grad0(
    const int num_of_quads,
    const scalar_t d42,
    const __global scalar2_t *pnts,
    const __global int4 *quads,
    global scalar2_t *df_aux)
{
    int tid = get_global_id(0);
    int id = tid*4;
    
    scalar2_t p0,p1,p2,p3;
    scalar2_t grad0, grad1, grad2, grad3;
    int4 f;

    if( tid < num_of_quads ) 
    {
        int4 f = quads[tid];
        p0 = pnts[ f.s0 ]; 
        p1 = pnts[ f.s1 ]; 
        p2 = pnts[ f.s2 ]; 
        p3 = pnts[ f.s3 ]; 

        grad0 = (scalar2_t)(0.0,0.0);
        grad1 = (scalar2_t)(0.0,0.0);
        grad2 = (scalar2_t)(0.0,0.0);
        grad3 = (scalar2_t)(0.0,0.0);
        
        imr_right_triangle_grad(p0,p1,p3,&grad0, &grad1, &grad3, d42);
        imr_right_triangle_grad(p1,p2,p0,&grad1, &grad2, &grad0, d42);
        imr_right_triangle_grad(p2,p3,p1,&grad2, &grad3, &grad1, d42);
        imr_right_triangle_grad(p3,p0,p2,&grad3, &grad0, &grad2, d42);
        
        df_aux[ id     ] = grad0;
        df_aux[ id + 1 ] = grad1;
        df_aux[ id + 2 ] = grad2;
        df_aux[ id + 3 ] = grad3;
    }
}

kernel void grad1(
         const int num_of_inner_verts,
         const global int *ia,
         const global int *ja,
         const global scalar2_t *df_aux,
         const scalar_t scale_factor,
         global scalar2_t *df)
{
    int tid = get_global_id(0);
    if( tid < num_of_inner_verts )
    {
        int j = tid;
        int end = ia[tid];

        scalar2_t d = (scalar2_t)(0.0,0.0);
        while( j < end ) {
            d += df_aux[ ja[j] ];
            j += num_of_inner_verts;
        }
        df[tid] = d * scale_factor;
    }
}

#if 0
kernel void grad1_laplace(
         const int num_of_inner_verts,
         const global int *ia,
         const global int *ja,
         const global scalar2_t *df_aux,
         const scalar_t scale_factor,
         global scalar2_t *df)
{
    int tid = get_global_id(0);
    if( tid < num_of_inner_verts )
    {
        int j = tid;
        int end = ia[tid];

        scalar2_t d = (scalar2_t)(0.0,0.0);
        int count = 0.0;
        while( j < end ) {
            d += df_aux[ 2*ja[j] ];
            j += num_of_inner_verts;
            ++count;
        }
        df[tid] = d / ((1.0 ) * count);
    }
}
#endif

scalar_t imr_triangle(
       const scalar2_t p0,
       const scalar2_t p1, 
       const scalar2_t p2,
       const scalar_t d42)
{
    scalar4_t a;    
    scalar_t h;

    a.even = p1 - p0;
    a.odd  =(2.0*p2 - p1 - p0) * M_SQRT1_3;
    h = (a.s0 * a.s3) - (a.s1 * a.s2);
    return dot(a,a) / (h + sqrt((h*h) + d42));
}


//!MIPS computation without the delta
scalar_t mips_triangle(
       const scalar2_t p0,
       const scalar2_t p1, 
       const scalar2_t p2,
       const scalar_t cota,
       const scalar_t cotb)
{
    scalar4_t a;
    scalar_t h;
    a.even = p1 - p0;
    a.odd = (p2 - p1)*cota + (p2 - p0)*cotb;
    h = (a.s0 * a.s3) - (a.s1 * a.s2);
    h = max(1e-6,h);  
    return dot(a,a) / (2.0 * h);
}

//!MIPS computation with delta
scalar_t mips_triangle_d(
       const scalar2_t p0,
       const scalar2_t p1, 
       const scalar2_t p2,
       const scalar_t cota,
       const scalar_t cotb,
       const scalar_t d42)
{
    scalar4_t a;
    scalar_t h;
    a.even = p1 - p0;
    a.odd = (p2 - p1)*cota + (p2 - p0)*cotb;
    //a.even = p1 - p0 - cota*(p2-p0);
    //a.odd = (p2 - p0)*(cota + cotb);
    h = (a.s0 * a.s3) - (a.s1 * a.s2);
    return dot(a,a) / (h + sqrt((h*h) + d42));
}

scalar_t imr_right_triangle(
       const scalar2_t p0,
       const scalar2_t p1, 
       const scalar2_t p2,
       const scalar_t d42)
{
    //a00, a10, a01, a11
    scalar4_t a;    

    a.lo = p1.xy - p0.xy;
    a.hi = p2.xy - p0.xy;
   
    scalar_t frobenius = dot(a,a);
    scalar_t h = (a.s0 * a.s3) - (a.s2 * a.s1);
    h = h + sqrt(h*h + d42);
    return frobenius / h;
}

kernel void compute_delta(
    const int num_of_tris,
    const int num_of_quads,
    const global scalar2_t *pnts,
    const global int4 *tris,
    const global int4 *quads,
    local scalar_t *work,
    global scalar_t *delta)
{
    int tid = get_local_id(0);
    int lsize = get_local_size(0);
    scalar4_t a;    
    int4 f;
    scalar2_t p0, p1, p2, p3;
    scalar_t min_s = 0.0;

    for( int i = tid; i < num_of_tris;  i += lsize )
    {
        f = tris[i];
        p0 = pnts[ f.s0 ]; 
        p1 = pnts[ f.s1 ]; 
        p2 = pnts[ f.s2 ]; 

        a.even = p1 - p0;
        a.odd  = (2.0*p2 - p1 - p0) * M_SQRT1_3;
        min_s = min( (a.s0 * a.s3) - (a.s2 * a.s1), min_s );
    }

    for( int i = tid; i < num_of_quads; i += lsize ) 
    {
        f = quads[i];
        p0 = pnts[ f.s0 ]; 
        p1 = pnts[ f.s1 ]; 
        p2 = pnts[ f.s2 ]; 
        p3 = pnts[ f.s3 ]; 
    
        a.lo = p1 - p0; a.hi = p3 - p0;
        min_s = min( (a.s0 * a.s3) - (a.s2 * a.s1), min_s );
        a.lo = p2 - p1; a.hi = p0 - p1;
        min_s = min( (a.s0 * a.s3) - (a.s2 * a.s1), min_s );
        a.lo = p3 - p2; a.hi = p1 - p2;
        min_s = min( (a.s0 * a.s3) - (a.s2 * a.s1), min_s );
        a.lo = p0 - p3; a.hi = p2 - p3;
        min_s = min( (a.s0 * a.s3) - (a.s2 * a.s1), min_s );
    }

    work[ tid ] = min_s;
    barrier(CLK_LOCAL_MEM_FENCE);

    //do reduction in shared mem
    for( int s = lsize/2; s > 0; s/=2  ) 
    {
        if( tid < s && (work[tid + s] < work[tid]) )
            work[tid] = work[tid + s];
        
        barrier(CLK_LOCAL_MEM_FENCE);
    }

    if(tid == 0)
        //delta[0] = (work[0] >= 0.0f ? 0.0f : sqrt(1e-04f*(1e-04f - work[0])));     
        //delta[0] = (work[0] >= 0.0 ? 0.0f : sqrt(1e-4*(1e-04 - work[0])));     
        //delta[0] = (work[0] >= 0.0 ? 0.0f : sqrt(1e-3*(1e-3 - work[0])));     
        //delta[0] = (work[0] >= 0.0 ? 0.0f : sqrt(1e-4*(1e-4 - work[0])));     
        delta[0] = (work[0] >= 0.0 ? 0.0f : sqrt(1e-5*(1e-5 - work[0])));     
        //delta[0] = (work[0] >= 0.0 ? 0.0f : sqrt(1e-3*(1e-3 - work[0])));     
}

kernel void line_search_prepare(
    const int n,
    const scalar_t alpha,
    const global scalar2_t *pnts,
    const global scalar2_t *df,
    global scalar2_t *pnts_aux)
{
    int tid = get_global_id(0);
    if( tid < n )
        pnts_aux[tid] = mad( alpha, df[tid], pnts[tid] );
}

kernel void mips_tris_f_d(
    const int num_of_triangles,
    const scalar_t d42,
    const global scalar2_t *pnts,
    const global int4 *triangles,
    const global scalar4_t *angles,
    global scalar_t *result)
{
    int tid = get_global_id(0);
    scalar2_t p0,p1,p2;
    int4 f;

    //compute knupp sum of all the triangles
    if( tid < num_of_triangles )
    {
        f = triangles[tid];
        p0 = pnts[ f.s0 ]; 
        p1 = pnts[ f.s1 ]; 
        p2 = pnts[ f.s2 ]; 
        
        result[tid] = mips_triangle_d(p0,p1,p2,angles[tid].s0,angles[tid].s1,d42);
    }
}

kernel void mips_tris_f(
    const int num_of_triangles,
    const scalar_t d42,
    const global scalar2_t *pnts,
    const global int4 *triangles,
    const global scalar4_t *angles,
    global scalar_t *result)
{
    int tid = get_global_id(0);
    scalar2_t p0,p1,p2;

    //compute knupp sum of all the triangles
    if( tid < num_of_triangles )
    {
        int4 f = triangles[tid];
        scalar2_t a = angles[tid].lo;

        p0 = pnts[ f.s0 ]; 
        p1 = pnts[ f.s1 ]; 
        p2 = pnts[ f.s2 ]; 
        
        result[tid] = mips_triangle(p0,p1,p2,a.s0,a.s1);
    }
}

kernel void combined_factor(
    const int num_of_triangles,
    const global scalar2_t *pnts,
    const global int4 *triangles,
    const global scalar4_t *angles,
    global scalar_t *result)
{
    int tid = get_global_id(0);

    //compute knupp sum of all the triangles
    if( tid < num_of_triangles )
    {
        int4 f = triangles[tid];
        scalar4_t a = angles[tid];
				scalar_t  area;

        scalar2_t p0 = pnts[ f.s0 ]; 
        scalar2_t p1 = pnts[ f.s1 ]; 
        scalar2_t p2 = pnts[ f.s2 ]; 

        scalar2_t e0 = p1 - p0;
        scalar2_t e1 = p2 - p1;
        scalar2_t e2 = p2 - p0;

        scalar_t det_w = max( e0.s0 * e2.s1 - e0.s1 * e2.s0, 1e-12 );
        result[tid] = sqrt( a.s3 / det_w ); 
    }
}

kernel void untangle_un(
    const int num_of_triangles,
    const global float2 *pnts,
    const global int4 *triangles,
    global float2 *df_aux,
    global float *areas)
{
    int tid = get_global_id(0);
    int id = tid*3;

    if( tid < num_of_triangles ) 
    {
        int4 f = triangles[tid];
        float2 p0 = pnts[ f.s0 ]; 
        float2 p1 = pnts[ f.s1 ]; 
        float2 p2 = pnts[ f.s2 ]; 

        float2 e0 = p1 - p0;
        float2 e1 = p2 - p1;
        float2 e2 = p0 - p2;
        double  det_w = -e0.s0 * e2.s1 + e0.s1 * e2.s0;
        e1 = (float2)(-e1.s1,e1.s0);
        e2 = (float2)(-e2.s1,e2.s0);
        e0 = (float2)(-e0.s1,e0.s0);
       
        if( det_w <= 0.0 ) {
          df_aux[ id     ] = (p1 - p0) + e1;
          df_aux[ id + 1 ] = (p2 - p1) + e2;
          df_aux[ id + 2 ] = (p0 - p2) + e0;
          //df_aux[ id     ] = +e1;
          //df_aux[ id + 1 ] = +e2;
          //df_aux[ id + 2 ] = +e0;
          det_w = -det_w;
        } else {
          df_aux[ id     ] = (float2)(0.0);
          df_aux[ id + 1 ] = (float2)(0.0);
          df_aux[ id + 2 ] = (float2)(0.0);
        }

#if 0
        if( det_w <= 0.0 ) {
          df_aux[ id     ] = e1;
          df_aux[ id + 1 ] = e2;
          df_aux[ id + 2 ] = e0;
          det_w = -det_w;
        } else {
          //df_aux[ id     ] = (float)(-2.0 * (det_w - 1e-10)) * e1;
          //df_aux[ id + 1 ] = (float)(-2.0 * (det_w - 1e-10)) * e2;
          //df_aux[ id + 2 ] = (float)(-2.0 * (det_w - 1e-10)) * e0;
          df_aux[ id     ] = (float2)(0.0);
          df_aux[ id + 1 ] = (float2)(0.0);
          df_aux[ id + 2 ] = (float2)(0.0);
        }
#endif

        areas[tid] = (float)fabs(det_w);
    }
}

kernel void grad1_un(
         const int num_of_inner_verts,
         global float2 *pnts,
         const global int *ia,
         const global int *ja,
         const global float2 *df_aux,
         const global float *areas,
         const float factor)
{
    int tid = get_global_id(0);
    if( tid < num_of_inner_verts )
    {
        int j = tid;
        int end = ia[tid];
        float2 untangle = (float2)(0.0,0.0);
        int count = 0;

        while( j < end ) {
            untangle += df_aux[ ja[j] ];
            j += num_of_inner_verts;
            ++count;
        }

        untangle /= ((float)count);
        pnts[tid] += factor * untangle;
    }
}

kernel void grad1_laplace(
         const int num_of_inner_verts,
         global float2 *pnts,
         const global int *ia,
         const global int *ja,
         const global float2 *df_aux,
         const global float *areas,
         const float factor)
{
    int tid = get_global_id(0);
    if( tid < num_of_inner_verts )
    {
        int j = tid;
        int end = ia[tid];
        float2 p = pnts[ tid ];
        float2 smooth = p * 1e-06f;
        float sumw    = 1e-06;

        while( j < end ) {
            int id = ja[j];
            float w = areas[ id / 3 ];
            sumw   += w; 
            smooth += w * df_aux[ id ];
            j += num_of_inner_verts;
        }

        pnts[tid] = mix( p, smooth / sumw, factor);
    }
}

kernel void untangle_laplace(
    const int num_of_triangles,
    const global float2 *pnts,
    const global int4 *triangles,
    global float2 *df_aux,
    global float *areas)
{
    int tid = get_global_id(0);
    int id = tid*3;

    if( tid < num_of_triangles ) 
    {
        int4 f = triangles[tid];

        float2 p0 = pnts[ f.s0 ]; 
        float2 p1 = pnts[ f.s1 ]; 
        float2 p2 = pnts[ f.s2 ]; 

        float2 e0 = p1 - p0;
        float2 e2 = p0 - p2;
        float det_w = e0.s0 * e2.s1 - e0.s1 * e2.s0;
        areas[tid] = fabs( det_w );
        df_aux[ id     ] = p1;
        df_aux[ id + 1 ] = p2;
        df_aux[ id + 2 ] = p0;
    }
}

kernel void combined_tris_untangle_grad0(
    const int num_of_triangles,
    const scalar_t d42,
    const global scalar2_t *pnts,
    const global int4 *triangles,
    const global scalar4_t *angles,
    global scalar2_t *df_aux)
{
    int tid = get_global_id(0);
    int id = tid*3;

    if( tid < num_of_triangles ) 
    {
        int4 f = triangles[tid];
        scalar4_t a = angles[tid]; // a.s3 is the area of original triangle
        a.s3 = a.s3 * a.s3;

        scalar2_t grad_0, grad_1, grad_2;

        scalar2_t p0 = pnts[ f.s0 ]; 
        scalar2_t p1 = pnts[ f.s1 ]; 
        scalar2_t p2 = pnts[ f.s2 ]; 
        scalar2_t e0 = p1 - p0;
        scalar2_t e1 = p2 - p1;
        scalar2_t e2 = p2 - p0;
#if UNTANGLE_FUNC1 

        scalar_t det_w = e0.s0 * e2.s1 - e0.s1 * e2.s0;

        scalar_t x0 = p0.s0;
        scalar_t y0 = p0.s1;
        scalar_t x1 = p1.s0;
        scalar_t y1 = p1.s1;
        scalar_t x2 = p2.s0;
        scalar_t y2 = p2.s1;

        //F = (Abs(det_w - b) - (det_w - b)) ^ 2
        if( det_w <= b ) {
          grad_0.s0 = -4.*( y1 - y2)*(2*PBETA - 2.*((x0 - x2)*(-y0 + y1) + (-x0 + x1)*(-y0 + y2)));
          grad_0.s1 = -4.*(-x1 + x2)*(2*PBETA - 2.*((x0 - x2)*(-y0 + y1) + (-x0 + x1)*(-y0 + y2)));
          grad_1.s0 = -4.*(-y0 + y2)*(2*PBETA - 2.*((x0 - x2)*(-y0 + y1) + (-x0 + x1)*(-y0 + y2)));
          grad_1.s1 = -4.*( x0 - x2)*(2*PBETA - 2.*((x0 - x2)*(-y0 + y1) + (-x0 + x1)*(-y0 + y2)));
          grad_2.s0 = -4.*( y0 - y1)*(2*PBETA - 2.*((x0 - x2)*(-y0 + y1) + (-x0 + x1)*(-y0 + y2)));
          grad_2.s1 = -4.*(-x0 + x1)*(2*PBETA - 2.*((x0 - x2)*(-y0 + y1) + (-x0 + x1)*(-y0 + y2)));
        } else {
          grad_0 = (scalar2_t)(0.0,0.0);
          grad_1 = (scalar2_t)(0.0,0.0);
          grad_2 = (scalar2_t)(0.0,0.0);
        }

	      df_aux[ id     ] = grad_0 / a.s3;
        df_aux[ id + 1 ] = grad_1 / a.s3;
        df_aux[ id + 2 ] = grad_2 / a.s3;
#else 
        scalar_t dir = a.s2*dot(e0,e0) + a.s1*dot(e2,e2) + a.s0*dot(e1,e1);
        scalar_t det_w = e0.s0 * e2.s1 - e0.s1 * e2.s0;
        scalar_t un_energy = ( (det_w - PBETA) * (det_w - PBETA) ) / a.s3;

        //F = (det_w - b) ^ 2
        scalar_t factor = d42 * ((2.0*(-PBETA -e2.s0*e0.s1 + e0.s0*e2.s1)) / a.s3);
        //grad_0.s0 = -e1.s1*factor;
        //grad_0.s1 =  e1.s0*factor;
        //grad_1.s0 =  e2.s1*factor;
        //grad_1.s1 = -e2.s0*factor;
        //grad_2.s0 = -e0.s1*factor;
        
        //grad_2.s1 =  e0.s0*factor;
        grad_0 = -2.0*((e0*a.s2 + e2*a.s1));
        grad_1 = 2.0*((e0*a.s2  - e1*a.s0));	      
        grad_2 = 2.0*((e1*a.s0  + e2*a.s1));        

        e1 = (scalar2_t)(-e1.s1,e1.s0);
        e2 = (scalar2_t)(e2.s1,-e2.s0);
        e0 = (scalar2_t)(-e0.s1,e0.s0);

        if( det_w < UN_LIMIT ) {
          df_aux[ id     ] = DIR_WEIGHT * grad_0 + (e1 * factor);
          df_aux[ id + 1 ] = DIR_WEIGHT * grad_1 + (e2 * factor);
          df_aux[ id + 2 ] = DIR_WEIGHT * grad_2 + (e0 * factor);
        } else {
          df_aux[ id     ] = DIR_WEIGHT * grad_0;
          df_aux[ id + 1 ] = DIR_WEIGHT * grad_1;
          df_aux[ id + 2 ] = DIR_WEIGHT * grad_2; 
        }
        
        //grad_0.s0 = 2*( y1 - y2)*(-b + (x0 - x2)*(-y0 + y1) + (-x0 + x1)*(-y0 + y2));
        //grad_0.s1 = 2*(-x1 + x2)*(-b + (x0 - x2)*(-y0 + y1) + (-x0 + x1)*(-y0 + y2));
        //grad_1.s0 = 2*(-y0 + y2)*(-b + (x0 - x2)*(-y0 + y1) + (-x0 + x1)*(-y0 + y2));
        //grad_1.s1 = 2*( x0 - x2)*(-b + (x0 - x2)*(-y0 + y1) + (-x0 + x1)*(-y0 + y2));
        //grad_2.s0 = 2*( y0 - y1)*(-b + (x0 - x2)*(-y0 + y1) + (-x0 + x1)*(-y0 + y2));
        //grad_2.s1 = 2*(-x0 + x1)*(-b + (x0 - x2)*(-y0 + y1) + (-x0 + x1)*(-y0 + y2));
#endif
    }
}


kernel void combined_tris_untangle_f(
    const int num_of_triangles,
    const scalar_t d42,
    const global scalar2_t *pnts,
    const global int4 *triangles,
    const global scalar4_t *angles,
    global scalar_t *result)
{
  int tid = get_global_id(0);

  if( tid < num_of_triangles ) {
        int4 f = triangles[tid];
        scalar4_t a = angles[tid];
        a.s3 = a.s3 * a.s3;

        scalar2_t p0 = pnts[ f.s0 ]; 
        scalar2_t p1 = pnts[ f.s1 ]; 
        scalar2_t p2 = pnts[ f.s2 ]; 

        scalar2_t e0 = p1 - p0;
        scalar2_t e1 = p2 - p1;
        scalar2_t e2 = p2 - p0;
        //scalar_t dir = a.s2*dot(e0,e0) + a.s1*dot(e2,e2) + a.s0*dot(e1,e1);
        scalar_t det_w = e0.s0 * e2.s1 - e0.s1 * e2.s0;
#if UNTANGLE_FUNC1
        result[tid] = pow( fabs(det_w - PBETA) - (det_w - PBETA), 2.0);
#else
        scalar_t dir = a.s2*dot(e0,e0) + a.s1*dot(e2,e2) + a.s0*dot(e1,e1);

        //F = (det_w - b) ^ 2
        if( det_w < UN_LIMIT ) {
          result[tid] =  DIR_WEIGHT * dir + d42 * (( (det_w - PBETA) * (det_w - PBETA) ) / a.s3);
        } else {
          result[tid] = DIR_WEIGHT * dir;
        }
         
#endif   
  }
}


kernel void combined_tris_f(
    const int num_of_triangles,
    const scalar_t d42,
    const global scalar2_t *pnts,
    const global int4 *triangles,
    const global scalar4_t *angles,
    global scalar_t *result)
{
    int tid = get_global_id(0);

    //compute knupp sum of all the triangles
    if( tid < num_of_triangles )
    {
        int4 f = triangles[tid];
        scalar4_t a = angles[tid];

        scalar2_t p0 = pnts[ f.s0 ]; 
        scalar2_t p1 = pnts[ f.s1 ]; 
        scalar2_t p2 = pnts[ f.s2 ]; 

#if 1
        scalar2_t e0 = p1 - p0;
        scalar2_t e1 = p2 - p1;
        scalar2_t e2 = p2 - p0;
        scalar_t det_w = max( e0.s0 * e2.s1 - e0.s1 * e2.s0, 1e-14 );
        scalar_t dir = a.s2*dot(e0,e0) + a.s1*dot(e2,e2) + a.s0*dot(e1,e1);
        result[tid] = dir / a.s3 + (dir * a.s3 / (det_w * det_w));
#else
        //scalar4_t ma;
        //ma.even = p1 - p0;
        //ma.odd = (p2 - p1)*a.s0 + (p2 - p0)*a.s1;
        //scalar_t det_w = (ma.s0 * ma.s3) - (ma.s1 * ma.s2);
        //det_w = max( det_w, 1e-8*a.s3 );
        //scalar_t dir = dot(ma,ma);
        //scalar4_t ma;
        scalar2_t e0 = p1 - p0;
        scalar2_t e1 = (p2 - p1)*a.s0 + (p2 - p0)*a.s1;
        scalar_t det_w = (e0.s0 * e1.s1) - (e1.s0 * e0.s1);
        det_w = max( det_w, a.s3*1e-10 );
        scalar_t dir = dot(e0,e0) + dot(e1,e1);
        //scalar_t det_w = (ma.s0 * ma.s3) - (ma.s1 * ma.s2);
        //det_w = max( det_w, 1e-8*a.s3 );
        //scalar_t dir = dot(ma,ma);
        
        result[tid] = dir * a.s2 + (dir / ((a.s2 * det_w) * det_w));
#endif

        //result[tid] = 0.5 * (dir/a.s3 + (dir*a.s3/(area*area)));

        //result[tid] = 0.5 * (dir/a.s3 + ((dir/area) * (a.s3/area)));

        //result[tid] = (( a.s2*(p1p0.s0*p1p0.s0 + p1p0.s1*p1p0.s1) +
        //                 a.s1*(p2p0.s0*p2p0.s0 + p2p0.s1*p2p0.s1) +
        //                 a.s0*(p2p1.s0*p2p1.s0 + p2p1.s1*p2p1.s1))/(2.*area)) *
        //                 a.s3; 

				//area = (p1.s0 - p0.s0)*(p2.s1 - p0.s1) - (p1.s1 - p0.s1)*(p2.s0 - p0.s0);
        //result[tid] = (( a.s2*((p1.s0 - p0.s0)*(p1.s0 - p0.s0) + (p1.s1 - p0.s1)*(p1.s1 - p0.s1)) +
        //                 a.s1*((p0.s0 - p2.s0)*(p0.s0 - p2.s0) + (p0.s1 - p2.s1)*(p0.s1 - p2.s1)) +
        //                 a.s0*((p2.s0 - p1.s0)*(p2.s0 - p1.s0) + (p2.s1 - p1.s1)*(p2.s1 - p1.s1)) )/(2.*area)) *
        //                 powr( (a.s3/area) + (area/a.s3), PARAMETERIZATION_THETA); 
    }
}

kernel void tris_f(
    const int num_of_triangles,
    const scalar_t d42,
    const global scalar2_t *pnts,
    const global int4 *triangles,
    global scalar_t *result)
{
    int tid = get_global_id(0);
    scalar2_t p0,p1,p2;
    int4 f;

    //compute knupp sum of all the triangles
    if( tid < num_of_triangles )
    {
        f = triangles[tid];
        p0 = pnts[ f.s0 ]; 
        p1 = pnts[ f.s1 ]; 
        p2 = pnts[ f.s2 ]; 
        
        result[tid] = imr_triangle(p0,p1,p2,d42);
    }
}

kernel void quads_f(
    const int num_of_quads,
    const scalar_t d42,
    const global scalar2_t *pnts,
    const global int4 *quads,
    global scalar_t *result)
{
    int tid = get_global_id(0);
    scalar2_t p0,p1,p2,p3;
    scalar4_t q;
    int4 f;

    //compute knupp sum of all the triangles
    if( tid < num_of_quads )
    {
        f = quads[tid];
        p0 = pnts[ f.s0 ]; 
        p1 = pnts[ f.s1 ]; 
        p2 = pnts[ f.s2 ]; 
        p3 = pnts[ f.s3 ]; 
        
        q.s0 = imr_right_triangle(p0,p1,p3,d42);
        q.s1 = imr_right_triangle(p1,p2,p0,d42);
        q.s2 = imr_right_triangle(p2,p3,p1,d42);
        q.s3 = imr_right_triangle(p3,p0,p2,d42);
        result[tid] = dot(q,q);
    }
}


