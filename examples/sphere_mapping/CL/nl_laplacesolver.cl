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

__constant float omega = 1.0;

__kernel void solve_conformal( 
        const __global float4 *v,
        const __global float4 *vf,
        const __global int *bounds,
        const __global int *verts,
        __global float4 *v_s,      //output
        const int numOfVertices,
        const __global float *w
        )
{
    //Vertex ID
    int id = get_global_id(0);

    if( id >= numOfVertices )
        return;

    float4 p = vf[id];   //vertex coordinates
    float4 c = (float4)(0.0f,0.0f,0.0f,0.0f);

    int start = bounds[id];
    int end = bounds[id+1];
    for( int i = start; i < end; ++i )
    {
        //neighbor vertex
        c += v[ verts[i] ] * w[i]; 
    }

    c = (1.0f - omega) * v[id] + omega*c;

    float l = dot(c,p) - 1.0f;
    v_s[id] = c - p*l;
}

__kernel void solve_conformal_res( 
        const __global float4 *v,
        const __global float4 *vf,
        const __global int *bounds,
        const __global int *verts,
        __global float4 *v_s,      //output
        const int numOfVertices,
        __global float *res,
        const __global float *w
        )
{
    //Vertex ID
    int id = get_global_id(0);

    if( id >= numOfVertices )
        return;

    float4 p = vf[id];   //vertex coordinates
    float4 c = (float4)(0.0f,0.0f,0.0f,0.0f);

    int start = bounds[id];
    int end = bounds[id+1];
    for( int i = start; i < end; ++i )
    {
        //neighbor vertex
        c += v[ verts[i] ] * w[i]; 
    }

    c = (1.0f - omega) * v[id] + omega*c;

    float l = dot(c,p) - 1.0f;
    float4 pos = c - p*l;
    v_s[id] = pos;
    res[id] = fast_distance(v[id],pos);
    //res[id] = dot(v[id]-pos,v[id]-pos);
}

__kernel void solve_equal( 
        const __global float4 *v,
        const __global float4 *vf,
        const __global int *bounds,
        const __global int *verts,
        __global float4 *v_s,      //output
        const int numOfVertices)
{
    //Vertex ID
    int id = get_global_id(0);

    if( id >= numOfVertices )
        return;

    float4 p = vf[id];   //vertex coordinates

    int start = bounds[id];
    int end = bounds[id+1];

    float4 c = v[verts[start]];
    for( int i = start+1; i < end; ++i )
        c += v[verts[i]]; 
    c /= (float)(end - start);
    c = (1.0f - omega) * v[id] + omega*c;

    float l = dot(c,p) - 1.0;
    v_s[id] = c - p*l;
}


__kernel void solve_equal_res( 
        const __global float4 *v,
        const __global float4 *vf,
        const __global int *bounds,
        const __global int *verts,
        __global float4 *v_s,      //output
        const int numOfVertices,
        __global float *res)
{
    //Vertex ID
    int id = get_global_id(0);

    if( id >= numOfVertices )
        return;

    float4 p = vf[id];   //vertex coordinates

    int start = bounds[id];
    int end = bounds[id+1];

    float4 c = v[verts[start]];
    for( int i = start+1; i < end; ++i )
        c += v[verts[i]]; 
    c /= (float)(end - start);
    c = (1.0f - omega) * v[id] + omega*c;

    float l = dot(c,p) - 1.0;
    float4 pos = c - p*l;
    v_s[id] = pos;
    res[id] = fast_distance(v[id],pos);
}

__kernel void normalize_solution(const __global float4 *v,  const int n, __global float4 *vf, __global float4 *v_s)
{
    int id = get_global_id(0);

    if( id >= n )
        return;

    vf[id] = v_s[id] = normalize(v[id]);
}

__kernel void convergence_conformal_res( 
        const __global float4 *v,
        const __global int *bounds,
        const __global int *verts,
        const int n,
        __global float *res,
        const __global float *w
        )
{
    //Vertex ID
    int id = get_global_id(0);

    if( id >= n )
        return;

    float4 c = (float4)(0.0f,0.0f,0.0f,0.0f);
 
    int start = bounds[id];
    int end = bounds[id+1];
    for( int i = start; i < end; ++i )
    {
        //neighbor vertex
        c += v[verts[i]] * w[i]; 
    }

    c = normalize(c);

    res[id] = dot(c - v[id],c - v[id]);
}

__kernel void convergence_equal_res( 
             const __global float4 *v,
             const __global int *bounds,
                     const __global int *verts,
                     const int n,
            __global float *res)        //residual output
{
    //Vertex ID
    int id = get_global_id(0);
        
    if( id >= n )
                return;
    
    int start = bounds[id];
    int end = bounds[id+1];
        
    float4 c = v[verts[start]];
        for( unsigned int i = start+1; i < end; ++i )
        c += v[verts[i]]; 
    c /= (float)(end - start);
    c = normalize(c);
    
    res[id] = dot(c - v[id],c - v[id]);
}

__kernel void lmax_residual(const int n, const const __global float *res, __global float *res_norma )
{
    float residue = 0.0;

    for( int i = 0; i < n; ++i )
        residue = max( res[i], residue );
   
    res_norma[0] = residue;
}

__kernel void l2_residual(const int n, const const __global float *res, __global float *res_norma )
{
    float residue = 0.0;

    for( int i = 0; i < n; ++i )
        residue += res[i];
   
    res_norma[0] = sqrt( residue ) / n;
}

