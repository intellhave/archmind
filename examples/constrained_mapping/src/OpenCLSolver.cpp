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

#include "OpenCLSolver.h"
#include "OpenCLBlas.h"
#include "NonLinearSolvers.h"
#include "UtilsCL.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <string>
#include <utility>
#include <cmath>
#include <boost/foreach.hpp>
#include <boost/bind.hpp>
#include <boost/function.hpp>
#include <boost/math/special_functions/fpclassify.hpp>

#ifndef isfinite
#define isfinite boost::math::isfinite
#endif

#ifndef isnan
#define isnan boost::math::isnan
#endif

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
using namespace parameterization;

// Round Up Division function
std::size_t round_up(int group_size, int global_size) 
{
    int r = global_size % group_size;
    if(r == 0) 
    {
        return global_size;
    } else 
    {
        return global_size + group_size - r;
    }
}

#define shft2(a,b,c) (a)=(b);(b)=(c);
#define shft3(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);
#define sign(a,b) ((b) >= 0.0 ? abs(a) : -abs(a))

typedef boost::function< cl_scalar_t (const cl_scalar_t &)> line_search_func_t;

template<typename Real,typename Func>
class base_line_search
{
public:
    base_line_search(Func &f, const Real lower_limit = -1.0, const Real upper_limit = 1e-13) : m_f(f)
    {
        m_LowerLimit = m_InitLower = lower_limit;
        m_AlphaHistory[0] = m_AlphaHistory[1] = m_AlphaHistory[2] = -1.0;
        m_UpperLimit = m_InitUpper = upper_limit;
        m_Fmin = 0.0;
    }

    Real line_search(Real &xmin, const Real &tol, const int max_iters = 50, const Real &wolf_acc = 0.0)
    {
        Real x0, x1, x2;
        Real fmin;
        x0 = m_LowerLimit;
        x1 = 0.0;

        //check Wolfe first condition
        fmin = m_f(-1.0);
        if( isfinite( fmin ) && (fmin <= m_Fmin - wolf_acc) ) {
            //std::cout << "wolfe_acc = " << wolf_acc << ", fmin = " << fmin << " prev = " << m_Fmin << "\n";
            m_Fmin = fmin;
            xmin = -1.0;
            return m_Fmin;
        }

        //if( fmin > m_Fmin - wolf_acc )
        //std::cout << "wolfe_acc = " << wolf_acc << ", fmin = " << fmin << " prev = " << m_Fmin << "\n";

        //while( !isfinite( m_f.eval_f(x0) ) ) {
        //    x0 *= 0.5;
        //    if( x0 > x1 - tol ) {
        //        m_Fmin = line_search_golden(x0,x1,x2,xmin,tol,max_iters);
        //        if( !isfinite(xmin) ){ xmin = 0.0; m_Fmin = 0.0;}
        //        return m_Fmin;
        //    }
        //}

        x2 = m_UpperLimit;
#if 0
        while( !isfinite( m_f.eval_f(x2) ) ) {
            x2 /= 10.0;
            if( x2 < x1 + tol ) {
                m_Fmin = line_search_golden(x0,x1,x2,xmin,tol,max_iters);
                if( !isfinite(xmin) ){ xmin = 0.0; m_Fmin = 0.0;}
                return m_Fmin;
            }
        }
#endif

        //std::cout << m_LowerLimit << "->" << x0 << "," << x1 << "," << m_UpperLimit << "->" << x2 << std::endl;
        m_Fmin = line_search_brent(x0,x1,x2,xmin,tol,max_iters);
        //m_Fmin = line_search_golden(x0,x1,x2,xmin,tol,max_iters);
        update_history(xmin);
        return m_Fmin;
    }

    Real line_search_brent(const Real &ax, const Real &bx, const Real &cx, Real &xmin, const Real &tol, const int max_iters = 50)
    {
        using std::abs;

        const Real zeps = 1.0e-12;
        const Real cgold = 0.3819660;
        Real a,b,d,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
        Real e=0.0;

        //std::cout << "ax = " << ax << " bx = " << bx << " cx = " << cx << "\n";
        //std::cout << "fax = " << eval_f(ax) << " fbx = " << eval_f(bx) << " fcx = " << eval_f(cx) << "\n";

        a = (ax < cx ? ax : cx);
        b = (ax > cx ? ax : cx);
        x = w = v = bx;

        //fw = fv = fx = m_f.eval_f(bx);
        if( m_Fmin == 0.0 ) m_Fmin = m_f(bx);
        fw = fv = fx = m_Fmin;

        for (int iter=1;iter<=max_iters;++iter) {
            xm=0.5*(a+b);
            tol2=2.0*(tol1=tol*abs(x)+zeps);
            if (abs(x-xm) <= (tol2-0.5*(b-a))) {
                xmin=x;
                return fx;
            }
            if (abs(e) > tol1) {
                r=(x-w)*(fx-fv);
                q=(x-v)*(fx-fw);
                p=(x-v)*q-(x-w)*r;
                q=2.0*(q-r);
                if (q > 0.0) p = -p;
                q=abs(q);
                etemp=e;
                e=d;
                if (abs(p) >= abs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x))
                    d=cgold*(e=(x >= xm ? a-x : b-x));
                else {
                    d=p/q;
                    u=x+d;
                    if (u-a < tol2 || b-u < tol2)
                        d=sign(tol1,xm-x);
                }
            } else {
                d=cgold*(e=(x >= xm ? a-x : b-x));
            }
            u=(abs(d) >= tol1 ? x+d : x+sign(tol1,d));
            fu=m_f(u);
            if (fu <= fx) {
                if (u >= x) a=x; else b=x;
                shft3(v,w,x,u)
                    shft3(fv,fw,fx,fu)
            } else {
                if (u < x) a=u; else b=u;
                if (fu <= fw || w == x) {
                    v=w;
                    w=u;
                    fv=fw;
                    fw=fu;
                } else if ( fu <= fv || v == x || v == w) {
                    v=u;
                    fv=fu;
                }
            }
        }
        xmin=x;
        return fx;
    }

    Real line_search_golden(const Real &ax, const Real &bx, const Real &cx, Real &xmin, const Real &tol, const int max_iters = 50)
    {
        using std::abs;

        const Real R = 0.61803399;
        const Real C = (1.0-R);
        Real f1,f2,x0,x1,x2,x3;

        x0=ax;
        x3=cx;
        if (abs(cx-bx) > abs(bx-ax)) {
            x1=bx;
            x2=bx+C*(cx-bx);
        } else {
            x2=bx;
            x1=bx-C*(bx-ax);
        }
        f1=m_f(x1);
        f2=m_f(x2);
        //std::cout << "f1 = " << f1 << "\n";
        //std::cout << "f2 = " << f2 << "\n";

        int iter = 0;
        while ((abs(x3-x0) > tol*(abs(x1)+abs(x2))) && ++iter < max_iters ) {
            if (isfinite(f2) && (f2 < f1)) {
                shft3(x0,x1,x2,R*x1+C*x3)
                    shft2(f1,f2,m_f(x2))
            } else {
                shft3(x3,x2,x1,R*x2+C*x0)
                shft2(f2,f1,m_f(x1))
            }
        }

        if (f1 < f2) {
            xmin=x1;
            return f1;
        } else {
            xmin=x2;
            return f2;
        }
    }

    void restart()
    {
        m_LowerLimit = m_InitLower;
        m_AlphaHistory[0] = m_AlphaHistory[1] = m_AlphaHistory[2] = -1.0;
        m_UpperLimit = m_InitUpper;
        m_Fmin = 0.0;
    }

private:
    void update_history(const Real &new_alpha)
    {
        using std::min;
        using std::max;

        shft3(m_AlphaHistory[0],m_AlphaHistory[1],m_AlphaHistory[2],new_alpha);
        Real avg = (m_AlphaHistory[0] + m_AlphaHistory[1] + m_AlphaHistory[2]) / 3.0;
        if( avg < 0.5 * m_LowerLimit ) 
          m_LowerLimit = max( -100, m_LowerLimit * 3.0 );
        else if( avg > 0.25 * m_LowerLimit )
          m_LowerLimit = min( m_LowerLimit / 2.0, -1e-06 );

        //std::cout << "(" << m_LowerLimit << "," << m_UpperLimit << ")\n";
    }

    Func &m_f;
    Real m_Fmin;
    Real m_AlphaHistory[3];
    Real m_LowerLimit;
    Real m_UpperLimit;
    Real m_InitLower;
    Real m_InitUpper;
};

template<typename Real>
double stable_sum( const int n, const Real *vec, bool &zero_area )
{
    Real result = 0.0;
    Real sum[4] = {0.0};
    Real c[4]   = {0.0};
    const int round4 = (n / 4) * 4;
    for( int i = 0; i < round4; i += 4 ) {
        Real y[4];
        Real t[4]; 

        if( vec[i + 0] == 0.0f ) zero_area = true;
        if( vec[i + 1] == 0.0f ) zero_area = true;
        if( vec[i + 2] == 0.0f ) zero_area = true;
        if( vec[i + 3] == 0.0f ) zero_area = true;

        y[0] = fabs(vec[i + 0]) - c[0];
        y[1] = fabs(vec[i + 1]) - c[1];
        y[2] = fabs(vec[i + 2]) - c[2];
        y[3] = fabs(vec[i + 3]) - c[3];

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

template<typename Real>
class mips_func
{
public:

    mips_func(cl::CommandQueue &q,
        int num_of_verts,
        int num_of_inner,
        int num_of_tris,
        int num_of_quads,
        cl::Buffer &pnts,
        cl::Buffer &pnts_aux,
        cl::Buffer &triangles,
        cl::Buffer &triangles_ia,
        cl::Buffer &triangles_ja,
        cl::Buffer &triangles_angles,
        cl::Buffer &triangles_factors,
        cl::Buffer &quads,
        cl::Buffer &quads_ia,
        cl::Buffer &quads_ja,
        cl::Buffer &df,
        double theta = 0.0) : 
    m_Queue(q),
        m_num_of_verts(num_of_verts),
        m_num_of_inner(num_of_inner),
        m_num_of_triangles(num_of_tris),
        m_num_of_quads(num_of_quads),
        m_Points(pnts),
        m_PointsAux(pnts_aux),
        m_Triangles(triangles),
        m_Triangles_ia(triangles_ia),
        m_Triangles_ja(triangles_ja),
        m_Triangle_angles(triangles_angles),
        m_Triangle_factors(triangles_factors),
        m_Quads(quads),
        m_Quads_ia(quads_ia),
        m_Quads_ja(quads_ja),
        m_D42(0.0),
        m_Df(df),
        m_blas_hd(q)
    {
        using namespace cl;
        using std::min;
        using std::max;

        //Read the source file
        ifstream source_file("constrained_solver.cl");
        if( !source_file ) 
            throw cl::Error(0, "Failed to open : \"constrained_solver.cl\"");

        try
        {
            m_Context = q.getInfo<CL_QUEUE_CONTEXT>();
            m_Device = q.getInfo<CL_QUEUE_DEVICE>();	

            m_num_of_threads = min( int(m_Device.getInfo<CL_DEVICE_MAX_WORK_GROUP_SIZE>()), 128 );
           
            std::string source_code(
                std::istreambuf_iterator<char>(source_file),
                (std::istreambuf_iterator<char>()));

            Program::Sources source(1, std::make_pair(source_code.c_str(), source_code.size()+1));
            Program program = Program(m_Context, source);

            std::vector< Device > devices;
            devices.push_back(m_Device);

            char options[255];
            try{
                //sprintf(options, "-DUSE_DOUBLE=%d -cl-mad-enable -cl-finite-math-only -cl-no-signed-zeros",CL_USE_DOUBLE);
                //sprintf(options, "-DUSE_DOUBLE=%d -cl-mad-enable -cl-finite-math-only -cl-no-signed-zeros -DPTHETA=%f -DTHETA_IS_ONE=%u",CL_USE_DOUBLE, theta, theta == 1.0 ? 1 : 0);
                sprintf(options, "-DUSE_DOUBLE=%d -cl-mad-enable -cl-no-signed-zeros -DPTHETA=%f -DTHETA_IS_ONE=%u",CL_USE_DOUBLE, theta, theta == 1.0 ? 1 : 0);
                std::cout << options << "\n";
                program.build(devices,options);
            } catch ( Error &err )
            {
                std::cerr << err.what() << "\n";
                std::cerr << program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(m_Device) << "\n";
                exit(1);
            }

            m_Result      = cl::Buffer(m_Context, CL_MEM_WRITE_ONLY, sizeof(cl_scalar_t) );
            m_TempBuffer  = cl::Buffer(m_Context, CL_MEM_READ_WRITE, sizeof(cl_scalar2_t) * max( 4*num_of_quads, 3*num_of_tris) );
            m_TempBuffer2 = cl::Buffer( m_Context, CL_MEM_READ_WRITE, sizeof(cl_scalar_t) * num_of_tris );
            m_TempBuffer3 = cl::Buffer( m_Context, CL_MEM_READ_WRITE, sizeof(cl_scalar2_t) * max( 4*num_of_quads, 3*num_of_tris)  );

            m_ComputeDelta = Kernel(program, "compute_delta");

            //compute_delta();
            set_delta(0.0);

            m_TotalScale = 1.0;
            if( theta <= 0.0 ) {
                std::cout << "Energy = Mips\n";
                m_EnergyType = 0;
                m_TrisFunc  = Kernel(program, "mips_tris_f");
                m_TrisGrad0 = Kernel(program, "mips_tris_grad0");
            } else {
                std::cout << "Energy = " << (m_D42 != 0.0 ? "Untangle\n" : "Isometric\n");
                m_EnergyType = m_D42 != 0.0 ? 2 : 1;
                m_TrisFunc  = Kernel(program, m_D42 != 0.0 ? "combined_tris_untangle_f"     : "combined_tris_f");
                m_TrisGrad0 = Kernel(program, m_D42 != 0.0 ? "combined_tris_untangle_grad0" : "combined_tris_grad0");
                m_D42 = 0.0;
            }

            m_QuadsFunc         = Kernel(program, "quads_f");
            m_LineSearchPrepare = Kernel(program, "line_search_prepare");
          
            m_QuadsGrad0        = Kernel(program, "quads_grad0");
            m_Grad1             = Kernel(program, "grad1");
            m_MipsFactor        = Kernel(program, "combined_factor");

            m_UnLaplaceIter     = Kernel(program, "untangle_laplace_d");
            m_UnLaplaceGrad1    = Kernel(program, "grad1_laplace_d");
            m_UnUntangleIter    = Kernel(program, "untangle_un_d");
            m_UnUntangleGrad1   = Kernel(program, "grad1_un_d");
        } catch ( Error &err )
        {
            std::cerr << "mips_func :: error in " << err.what() << " code = " << err.err() << "\n";
        }
    }

    double rescale()
    {
        if( m_EnergyType == 1 ) {
            Real scale = eval_factor();     //compute initial approximation
            Real initial = scale;
            line_search_func_t func = boost::bind( &mips_func<cl_scalar_t>::eval_scale_f, this, _1 );
            base_line_search<cl_scalar_t, line_search_func_t> kls( func );
            kls.line_search_brent(0.1, scale, 5.0*scale, scale, 1e-06, 30 ); 
            arch::blas::scal( m_blas_hd, m_num_of_verts*2, scale, m_Points, 1 );

            std::cout << "New Scale Factor = " << scale << " Initial = " << initial << "\n";
            m_TotalScale *= scale;
            return scale;
        } else if( m_EnergyType == 2 ) {
            if( m_D42 == 0.0 ) m_D42 = 0.1;
            else m_D42 *= 10.0;
            std::cout << "new delta = " << m_D42 << "\n";
        }

        return 1.0;
    }

    Real get_total_scale() {return m_TotalScale;}

    Real eval_factor()
    {
        double scale_factor = 1.0;
        try{
            cl::Event event;
            int reduction_threads = m_num_of_threads;

            if( m_num_of_triangles > 0 )
            {
                m_MipsFactor.setArg(0,m_num_of_triangles); 
                m_MipsFactor.setArg(1,m_Points);
                m_MipsFactor.setArg(2,m_Triangles);
                m_MipsFactor.setArg(3,m_Triangle_angles);
                m_MipsFactor.setArg(4,m_TempBuffer);

                m_Queue.enqueueNDRangeKernel(
                    m_MipsFactor,
                    cl::NullRange,
                    cl::NDRange( arch::clutils::round_up(m_num_of_threads, m_num_of_triangles) ),
                    cl::NDRange( m_num_of_threads),
                    NULL,
                    &event);

                event.wait();
               
                arch::blas::avg( m_blas_hd, m_num_of_triangles, m_TempBuffer, &scale_factor );
            }

        } catch ( cl::Error &err ) {
            std::cerr << "eval_factor :: error in " << err.what() << " code = " << err.err() << " error = " << arch::clutils::error(err.err() ) << "\n";
        }

        return scale_factor;
    }

    Real eval_scale_f(const Real &alpha)
    {
        ++m_Stats.fevals;
        scalar_t results_tris = 0.0;
        try{
            cl::Event event;
            arch::blas::copy( m_blas_hd, m_num_of_verts*2, m_Points, 1, m_PointsAux, 1);
            arch::blas::scal( m_blas_hd, m_num_of_verts*2, alpha, m_PointsAux, 1 );

            int reduction_threads = m_num_of_threads;

            if( m_num_of_triangles > 0 )
            {
                m_TrisFunc.setArg(0,m_num_of_triangles); 
                m_TrisFunc.setArg(1,m_D42); 
                m_TrisFunc.setArg(2,m_PointsAux);
                m_TrisFunc.setArg(3,m_Triangles);
                m_TrisFunc.setArg(4,m_Triangle_angles);
                //m_TrisFunc.setArg(4,m_Triangle_factors);
                m_TrisFunc.setArg(5,m_TempBuffer);

                m_Queue.enqueueNDRangeKernel(
                    m_TrisFunc,
                    cl::NullRange,
                    cl::NDRange( arch::clutils::round_up(m_num_of_threads, m_num_of_triangles) ),
                    cl::NDRange( m_num_of_threads),
                    NULL,
                    &event);

                event.wait();
               
                arch::blas::avg( m_blas_hd, m_num_of_triangles, m_TempBuffer, &results_tris );
                results_tris = results_tris - 1.0;
            }

        } catch ( cl::Error &err ) {
            std::cerr << "eval_scale_f :: error in " << err.what() << " code = " << err.err() << " error = " << arch::clutils::error(err.err() ) << "\n";
        }
        return results_tris;
    }

    Real eval_f(const Real &alpha = 0.0)
    {
        ++m_Stats.fevals;
        scalar_t result_quads = 0.0;
        scalar_t results_tris = 0.0;
        //std::cout << "Eval f\n";
        try{
            cl::Event event;

            if( alpha != 0.0 && alpha != -1.0 )
            {
                m_LineSearchPrepare.setArg(0,m_num_of_inner);
                m_LineSearchPrepare.setArg(1,alpha);
                m_LineSearchPrepare.setArg(2,m_Points);
                m_LineSearchPrepare.setArg(3,m_Df);
                m_LineSearchPrepare.setArg(4,m_PointsAux);

                m_Queue.enqueueNDRangeKernel(
                    m_LineSearchPrepare,
                    cl::NullRange,
                    cl::NDRange( arch::clutils::round_up(m_num_of_threads, m_num_of_inner) ),
                    cl::NDRange(m_num_of_threads),
                    NULL,
                    &event);

                if( PROFILE_OPENCL ) {
                    event.wait();
                    cl_long start_time = event.getProfilingInfo<CL_PROFILING_COMMAND_START>();
                    cl_long end_time = event.getProfilingInfo<CL_PROFILING_COMMAND_END>();
                    std::cout << "prepare time = " << ((end_time - start_time) / 1000.0) << "  microsecond\n";
                }
            } else if( alpha == -1.0 ) {
                arch::blas::xmy( m_blas_hd, m_num_of_inner*2, m_Points, m_Df, m_PointsAux );
            }

            int reduction_threads = m_num_of_threads;

            if( m_num_of_triangles > 0 )
            {
                m_TrisFunc.setArg(0,m_num_of_triangles); 
                m_TrisFunc.setArg(1,m_D42); 
                m_TrisFunc.setArg(2,alpha != 0.0 ? m_PointsAux : m_Points);
                m_TrisFunc.setArg(3,m_Triangles);
                m_TrisFunc.setArg(4,m_Triangle_angles);
                //m_TrisFunc.setArg(4,m_Triangle_factors);
                m_TrisFunc.setArg(5,m_TempBuffer);

                m_Queue.enqueueNDRangeKernel(
                    m_TrisFunc,
                    cl::NullRange,
                    cl::NDRange( arch::clutils::round_up(m_num_of_threads, m_num_of_triangles) ),
                    cl::NDRange(m_num_of_threads),
                    NULL,
                    &event);

                if( PROFILE_OPENCL ) {
                    event.wait();
                    cl_long start_time = event.getProfilingInfo<CL_PROFILING_COMMAND_START>();
                    cl_long end_time = event.getProfilingInfo<CL_PROFILING_COMMAND_END>();
                    std::cout << "eval_f time = " << ((end_time - start_time) / 1000.0) << "  microsecond\n";
                }
                  
                //arch::blas::sum( m_blas_hd, m_num_of_triangles, m_TempBuffer, m_TempBuffer, &results_tris );
                //results_tris = (results_tris / m_num_of_triangles) - 1.0;
                arch::blas::avg( m_blas_hd, m_num_of_triangles, m_TempBuffer, &results_tris );
                results_tris = results_tris - 1.0;

            }

            if( m_num_of_quads > 0 ) 
            {
                m_QuadsFunc.setArg(0,m_num_of_quads); 
                m_QuadsFunc.setArg(1,m_D42); 
                m_QuadsFunc.setArg(2,alpha != 0.0 ? m_PointsAux : m_Points);
                m_QuadsFunc.setArg(3,m_Quads);
                m_QuadsFunc.setArg(4,m_TempBuffer);

                m_Queue.enqueueNDRangeKernel(
                    m_QuadsFunc,
                    cl::NullRange,
                    cl::NDRange( arch::clutils::round_up(m_num_of_threads, m_num_of_quads) ),
                    cl::NDRange(m_num_of_threads));
                   
                //arch::blas::sum( m_blas_hd, m_num_of_quads, m_TempBuffer, m_TempBuffer, &result_quads );
                //result_quads = (result_quads / (4.0 * m_num_of_quads)) - 1.0;
                arch::blas::avg( m_blas_hd, m_num_of_quads, m_TempBuffer, &result_quads );
                result_quads = (result_quads / 4.0) - 1.0;
                //result_quads = (result_quads / 4.0);
            }
        } catch ( cl::Error &err ) {
            std::cerr << "eval_f :: error in " << err.what() << " code = " << err.err() << " error = " << arch::clutils::error(err.err() ) << "\n";
        }

        return results_tris + result_quads;
    }

    double eval_laplace_f(double smooth_factor, bool ret_value)
    {
        try{
            if( m_num_of_triangles > 0 ) 
            {
                m_UnLaplaceIter.setArg(0,m_num_of_triangles); 
                m_UnLaplaceIter.setArg(1,m_Points);
                m_UnLaplaceIter.setArg(2,m_Triangles);
                m_UnLaplaceIter.setArg(3,m_TempBuffer);
                m_UnLaplaceIter.setArg(4,m_TempBuffer2 );

                cl::Event event;
                m_Queue.enqueueNDRangeKernel(
                    m_UnLaplaceIter,
                    cl::NullRange,
                    cl::NDRange(round_up(m_num_of_threads,m_num_of_triangles) ),
                    cl::NDRange(m_num_of_threads),
                    NULL,
                    &event);  

                if( PROFILE_OPENCL ) {
                    event.wait();
                    unsigned start_time = event.getProfilingInfo<CL_PROFILING_COMMAND_START>();
                    unsigned end_time = event.getProfilingInfo<CL_PROFILING_COMMAND_END>();
                    std::cout << "kiter time = " << ((end_time - start_time) / 1000.0) << "  microsecond\n";
                }
                
                m_UnLaplaceGrad1.setArg(0,m_num_of_inner ); 
                m_UnLaplaceGrad1.setArg(1,m_Points ); 
                m_UnLaplaceGrad1.setArg(2,m_Triangles_ia);
                m_UnLaplaceGrad1.setArg(3,m_Triangles_ja);
                m_UnLaplaceGrad1.setArg(4,m_TempBuffer);
                m_UnLaplaceGrad1.setArg(5,m_TempBuffer2 );
                m_UnLaplaceGrad1.setArg(6,smooth_factor);
                
                m_Queue.enqueueNDRangeKernel(
                    m_UnLaplaceGrad1,
                    cl::NullRange,
                    cl::NDRange(round_up(64,m_num_of_inner) ),
                    cl::NDRange(64),
                    NULL,
                    &event);

                if( PROFILE_OPENCL ) {
                    event.wait();
                    unsigned start_time = event.getProfilingInfo<CL_PROFILING_COMMAND_START>();
                    unsigned end_time = event.getProfilingInfo<CL_PROFILING_COMMAND_END>();
                    std::cout << "kgrad time = " << ((end_time - start_time) / 1000.0) << "  microsecond\n";
                }

                if( ret_value ) {
                    std::vector< double > temp(m_num_of_triangles);
                    m_blas_hd.BlasQueue.enqueueReadBuffer( m_TempBuffer2, GL_TRUE, 0, sizeof(cl_double) * m_num_of_triangles, &temp[0] );
                    return arch::blas::stable_sum(m_num_of_triangles, &temp[0]) / 2.0;
                } else return 0.0;
            }
        } catch ( cl::Error &err ) {
            std::cerr << "eval_laplace_f :: error in " << err.what() << " code = " << err.err() << " error = " << arch::clutils::error(err.err() ) << "\n";
        }
    }

    double eval_untangle_f(double un_factor, bool ret_value)
    {
        try{
            if( m_num_of_triangles > 0 ) 
            {
                m_UnUntangleIter.setArg(0,m_num_of_triangles); 
                m_UnUntangleIter.setArg(1,m_Points);
                m_UnUntangleIter.setArg(2,m_Triangles);
                m_UnUntangleIter.setArg(3,m_TempBuffer);
                m_UnUntangleIter.setArg(4,m_TempBuffer3 );
                m_UnUntangleIter.setArg(5,m_TempBuffer2 );

                cl::Event event;
                m_Queue.enqueueNDRangeKernel(
                    m_UnUntangleIter,
                    cl::NullRange,
                    cl::NDRange(round_up(m_num_of_threads,m_num_of_triangles) ),
                    cl::NDRange(m_num_of_threads),
                    NULL,
                    &event);  

                if( PROFILE_OPENCL ) {
                    event.wait();
                    unsigned start_time = event.getProfilingInfo<CL_PROFILING_COMMAND_START>();
                    unsigned end_time = event.getProfilingInfo<CL_PROFILING_COMMAND_END>();
                    std::cout << "kiter time = " << ((end_time - start_time) / 1000.0) << "  microsecond\n";
                }
                

                //arch::blas::print( m_blas_hd, m_num_of_triangles * 3 * 2, m_TempBuffer );

                m_UnUntangleGrad1.setArg(0,m_num_of_inner ); 
                m_UnUntangleGrad1.setArg(1,m_Points ); 
                m_UnUntangleGrad1.setArg(2,m_Triangles_ia);
                m_UnUntangleGrad1.setArg(3,m_Triangles_ja);
                m_UnUntangleGrad1.setArg(4,m_TempBuffer);
                m_UnUntangleGrad1.setArg(5,m_TempBuffer3);
                m_UnUntangleGrad1.setArg(6,m_TempBuffer2 );
                m_UnUntangleGrad1.setArg(7,un_factor);
                
                m_Queue.enqueueNDRangeKernel(
                    m_UnUntangleGrad1,
                    cl::NullRange,
                    cl::NDRange(round_up(64,m_num_of_inner) ),
                    cl::NDRange(64),
                    NULL,
                    &event);

                if( PROFILE_OPENCL ) {
                    event.wait();
                    unsigned start_time = event.getProfilingInfo<CL_PROFILING_COMMAND_START>();
                    unsigned end_time = event.getProfilingInfo<CL_PROFILING_COMMAND_END>();
                    std::cout << "kgrad time = " << ((end_time - start_time) / 1000.0) << "  microsecond\n";
                }
                
                if( ret_value ) {
                    std::vector< double > temp(m_num_of_triangles);
                    m_blas_hd.BlasQueue.enqueueReadBuffer( m_TempBuffer2, CL_TRUE, 0, sizeof(cl_double) * m_num_of_triangles, &temp[0] );
                    bool zero_area;
                    double rev_val = stable_sum(m_num_of_triangles, &temp[0], zero_area) / 2.0;
                    //if( zero_area ) std::cout << "Zero area detected\n";
                    return rev_val;
                    //return arch::blas::stable_sum(m_num_of_triangles, &temp[0]) / 2.0;
                } 
            }
        } catch ( cl::Error &err ) {
            std::cerr << "eval_laplace_f :: error in " << err.what() << " code = " << err.err() << " error = " << arch::clutils::error(err.err() ) << "\n";
        }

        return 0.0;
    }

#if 0
    void eval_grad_f()
    {
        std::vector< cl_scalar_t > px(2*m_num_of_inner);
        std::vector< cl_scalar_t > df(2*m_num_of_inner);
        std::vector< cl_scalar_t > x(2*m_num_of_inner);

        m_Queue.enqueueReadBuffer(m_Points, CL_TRUE, 0, x.size()*sizeof(cl_scalar_t), &x[0] );

        for( int i = 0; i < 2*m_num_of_inner; ++i ) {
            for( int j = 0; j < 2*m_num_of_inner; ++j )
                px[j] = ((i == j) ? (x[j]-1e-8) : x[j]);
            m_Queue.enqueueWriteBuffer(m_Points, CL_TRUE, 0, px.size()*sizeof(cl_scalar_t), &px[0] );
            cl_scalar_t f0 = eval_f();
            for( int j = 0; j < 2*m_num_of_inner; ++j )
                px[j] = ((i == j) ? (x[j]+1e-8) : x[j]);
            m_Queue.enqueueWriteBuffer(m_Points, CL_TRUE, 0, px.size()*sizeof(cl_scalar_t), &px[0] );
            cl_scalar_t f1 = eval_f();

            //cl_scalar_t scale_factor = 1.0 / ( 4.0*m_num_of_quads + m_num_of_triangles);
            df[i] = (f1 - f0) / (2*1e-08);
        }

        m_Queue.enqueueWriteBuffer(m_Points, CL_TRUE, 0, x.size()*sizeof(cl_scalar_t), &x[0] );
        m_Queue.enqueueWriteBuffer(m_Df, CL_TRUE, 0, df.size()*sizeof(cl_scalar_t), &df[0] );
    }
#else
    void eval_grad_f()
    {
        ++m_Stats.gradevals; 
        //compute_delta();

        try{
            if( m_num_of_triangles > 0 ) 
            {
                m_TrisGrad0.setArg(0,m_num_of_triangles); 
                m_TrisGrad0.setArg(1,m_D42);
                m_TrisGrad0.setArg(2,m_Points);
                m_TrisGrad0.setArg(3,m_Triangles);
                m_TrisGrad0.setArg(4,m_Triangle_angles);
                m_TrisGrad0.setArg(5,m_TempBuffer);

                cl::Event event;
                m_Queue.enqueueNDRangeKernel(
                    m_TrisGrad0,
                    cl::NullRange,
                    cl::NDRange(round_up(m_num_of_threads,m_num_of_triangles) ),
                    cl::NDRange(m_num_of_threads),
                    NULL,
                    &event);

                if( PROFILE_OPENCL ) {
                    event.wait();
                    cl_long start_time = event.getProfilingInfo<CL_PROFILING_COMMAND_START>();
                    cl_long end_time = event.getProfilingInfo<CL_PROFILING_COMMAND_END>();
                    std::cout << "grad_f time = " << ((end_time - start_time) / 1000.0) << "  microsecond\n";
                }

                double scale_factor = 1.0 / m_num_of_triangles;
                m_Grad1.setArg(0,m_num_of_inner ); 
                m_Grad1.setArg(1,m_Triangles_ia);
                m_Grad1.setArg(2,m_Triangles_ja);
                m_Grad1.setArg(3,m_TempBuffer);
                m_Grad1.setArg(4,scale_factor);
                m_Grad1.setArg(5,m_Df);

                m_Queue.enqueueNDRangeKernel(
                    m_Grad1,
                    cl::NullRange,
                    cl::NDRange(round_up(64,m_num_of_inner) ),
                    cl::NDRange(64),
                    NULL,
                    &event);

                if( PROFILE_OPENCL ) {
                    event.wait();
                    cl_long start_time = event.getProfilingInfo<CL_PROFILING_COMMAND_START>();
                    cl_long end_time = event.getProfilingInfo<CL_PROFILING_COMMAND_END>();
                    std::cout << "grad_gather time = " << ((end_time - start_time) / 1000.0) << "  microsecond\n";
                }
            }

#if 0
            if( m_num_of_quads > 0 ) 
            {
                m_QuadsGrad0.setArg(0,m_num_of_quads); 
                m_QuadsGrad0.setArg(1,m_D42);
                m_QuadsGrad0.setArg(2,m_Points);
                m_QuadsGrad0.setArg(3,m_Quads);
                m_QuadsGrad0.setArg(4,m_TempBuffer);

                m_Queue.enqueueNDRangeKernel(
                    m_QuadsGrad0,
                    cl::NullRange,
                    cl::NDRange(round_up(m_num_of_threads,m_num_of_quads) ),
                    cl::NDRange(m_num_of_threads));

                m_Grad1.setArg(0,m_num_of_inner ); 
                m_Grad1.setArg(1,m_Quads_ia);
                m_Grad1.setArg(2,m_Quads_ja);
                m_Grad1.setArg(3,m_TempBuffer);
                m_Grad1.setArg(4,m_Df);

                m_Queue.enqueueNDRangeKernel(
                    m_Grad1,
                    cl::NullRange,
                    cl::NDRange(round_up(m_num_of_threads,m_num_of_inner) ),
                    cl::NDRange(m_num_of_threads));
            }
#endif
            //normalize the result
            //cl_scalar_t scale_factor = 1.0 / ( 4.0*m_num_of_quads + m_num_of_triangles);
            //arch::blas::scal( m_blas_hd, m_num_of_inner*2, scale_factor, m_Df, 1 );
        } catch ( cl::Error &err ) {
            std::cerr << "eval_grad_f :: error in " << err.what() << " code = " << err.err() << " error = " << arch::clutils::error(err.err() ) << "\n";
        }
    }
#endif

    void set_delta(const Real &delta)
    {
        m_D42 = 4.0 * delta * delta;
    }

    cl::Buffer &get_x() 
    {
        return m_Points;
    }

    cl::Buffer &get_df()
    {
        return m_Df;
    }

    struct stats
    {
        stats() : fevals(0), gradevals(0) {}
        cl_int fevals;
        cl_int gradevals;
    };

    stats get_stats()const { return m_Stats; }

    void check_derivatives(double eps = 1e-03)
    {
        std::cout << "Starting derivative checker...\n";
        std::vector< cl_scalar_t > px(2*m_num_of_inner);
        std::vector< cl_scalar_t > df(2*m_num_of_inner);
        std::vector< cl_scalar_t > x(2*m_num_of_inner);

        m_Queue.enqueueReadBuffer(m_Points, CL_TRUE, 0, x.size()*sizeof(cl_scalar_t), &x[0] );
        m_Queue.enqueueReadBuffer(m_Df, CL_TRUE, 0, df.size()*sizeof(cl_scalar_t), &df[0] );

        for( int i = 0; i < 2*m_num_of_inner; ++i ) {
            for( int j = 0; j < 2*m_num_of_inner; ++j )
                px[j] = ((i == j) ? (x[j]-1e-8) : x[j]);
            m_Queue.enqueueWriteBuffer(m_Points, CL_TRUE, 0, px.size()*sizeof(cl_scalar_t), &px[0] );
            cl_scalar_t f0 = eval_f();
            for( int j = 0; j < 2*m_num_of_inner; ++j )
                px[j] = ((i == j) ? (x[j]+1e-8) : x[j]);
            m_Queue.enqueueWriteBuffer(m_Points, CL_TRUE, 0, px.size()*sizeof(cl_scalar_t), &px[0] );
            cl_scalar_t f1 = eval_f();

            //cl_scalar_t scale_factor = 1.0 / ( 4.0*m_num_of_quads + m_num_of_triangles);
            cl_scalar_t fd = (f1 - f0) / (2*1e-08);

            if( (fabs( fd - df[i] ) > eps) ) 
            {
                printf("grad[%d] = %.12lf~%.12lf [ %.12lf ]\n", i, fd, df[i], fabs( fd - df[i] ) );
                //std::cout << "grad_f[" << i << "] = " << fd << "~" << df[i] << "\n";
            }
        }
        std::cout << "derivative checker finished\n";

        m_Queue.enqueueWriteBuffer(m_Points, CL_TRUE, 0, x.size()*sizeof(cl_scalar_t), &x[0] );
    }

    void compute_delta()
    {
        cl_scalar_t delta;

        try{
            int reduction_threads = m_num_of_threads;

            m_ComputeDelta.setArg(0,m_num_of_triangles); 
            m_ComputeDelta.setArg(1,m_num_of_quads);
            m_ComputeDelta.setArg(2,m_Points); 
            m_ComputeDelta.setArg(3,m_Triangles);
            m_ComputeDelta.setArg(4,m_Quads);
            m_ComputeDelta.setArg(5,sizeof(cl_scalar_t)*reduction_threads,NULL);
            m_ComputeDelta.setArg(6,m_Result);

            m_Queue.enqueueNDRangeKernel(
                m_ComputeDelta,
                cl::NullRange,
                cl::NDRange(reduction_threads),
                cl::NDRange(reduction_threads));

            m_Queue.enqueueReadBuffer(m_Result, CL_TRUE, 0, sizeof(scalar_t), &delta);
            if( delta != 0.0 ) std::cout << "Delta = " << delta << "\n";
            set_delta(delta);

        } catch ( cl::Error &err ) {
            std::cerr << "compute_delta :: error in " << err.what() << " code = " << err.err() << " error = " << arch::clutils::error(err.err() ) << "\n";
        }
    }

private:

    cl::Context m_Context;
    cl::Device m_Device;
    cl::CommandQueue &m_Queue;
    cl_int m_num_of_verts;
    cl_int m_num_of_inner;
    cl_int m_num_of_triangles;
    cl_int m_num_of_quads;
    cl_int m_num_of_threads;
    cl_scalar_t m_D42;
    cl::Buffer &m_Points;
    cl::Buffer &m_PointsAux;
    cl::Buffer &m_Triangles;
    cl::Buffer &m_Triangles_ia;
    cl::Buffer &m_Triangles_ja;
    cl::Buffer &m_Triangle_angles;
    cl::Buffer &m_Triangle_factors;
    cl::Buffer &m_Quads;
    cl::Buffer &m_Quads_ia;
    cl::Buffer &m_Quads_ja;
    cl::Buffer &m_Df;
    cl::Buffer m_Result;
    cl::Buffer m_TempBuffer;
    cl::Buffer m_TempBuffer2;
    cl::Buffer m_TempBuffer3;
    cl::Kernel m_TrisFunc;
    cl::Kernel m_QuadsFunc;
    cl::Kernel m_ComputeDelta;
    cl::Kernel m_LineSearchPrepare;
    cl::Kernel m_TrisGrad0;
    cl::Kernel m_QuadsGrad0;
    cl::Kernel m_Grad1;
    cl::Kernel m_MipsFactor;
    cl::Kernel m_UnLaplaceIter;
    cl::Kernel m_UnLaplaceGrad1;
    cl::Kernel m_UnUntangleIter;
    cl::Kernel m_UnUntangleGrad1;

    arch::blas::handler m_blas_hd;
    stats m_Stats;
    Real m_TotalScale;
    int m_EnergyType;
};

template<typename Real>
class knupp_func
{
public:

    knupp_func(cl::CommandQueue &q,
        int num_of_inner,
        int num_of_tris,
        int num_of_quads,
        cl::Buffer &pnts,
        cl::Buffer &pnts_aux,
        cl::Buffer &triangles,
        cl::Buffer &triangles_ia,
        cl::Buffer &triangles_ja,
        cl::Buffer &quads,
        cl::Buffer &quads_ia,
        cl::Buffer &quads_ja,
        cl::Buffer &df) : 
    m_Queue(q),
        m_num_of_inner(num_of_inner),
        m_num_of_triangles(num_of_tris),
        m_num_of_quads(num_of_quads),
        m_Points(pnts),
        m_PointsAux(pnts_aux),
        m_Triangles(triangles),
        m_Triangles_ia(triangles_ia),
        m_Triangles_ja(triangles_ja),
        m_Quads(quads),
        m_Quads_ia(quads_ia),
        m_Quads_ja(quads_ja),
        m_D42(0.0),
        m_Df(df),
        m_blas_hd(q)
    {
        using namespace cl;
        using std::min;
        using std::max;

         //Read the source file
         ifstream source_file("constrained_solver.cl");
         if( !source_file ) 
             throw cl::Error(0, "Failed to open : \"constrained_solver.cl\"");

        try
        {
            m_Context = q.getInfo<CL_QUEUE_CONTEXT>();
            m_Device = q.getInfo<CL_QUEUE_DEVICE>();	

            m_num_of_threads = min( int(m_Device.getInfo<CL_DEVICE_MAX_WORK_GROUP_SIZE>()), 64 );

            std::string source_code(
                std::istreambuf_iterator<char>(source_file),
                (std::istreambuf_iterator<char>()));

            Program::Sources source(1, std::make_pair(source_code.c_str(), source_code.size()+1));
            Program program = Program(m_Context, source);

            std::vector< Device > devices;
            devices.push_back(m_Device);

            char options[255];
            try{
                sprintf(options, "-DUSE_DOUBLE=%d  -cl-mad-enable -cl-fast-relaxed-math",CL_USE_DOUBLE);
                program.build(devices,options);
            } catch ( Error &err )
            {
                std::cerr << err.what() << "\n";
                std::cerr << program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(m_Device) << "\n";
                exit(1);
            }

            m_TrisFunc          = Kernel(program, "tris_f");
            m_QuadsFunc         = Kernel(program, "quads_f");
            m_LineSearchPrepare = Kernel(program, "line_search_prepare");
            m_ComputeDelta      = Kernel(program, "compute_delta");
            m_TrisGrad0         = Kernel(program, "tris_grad0");
            m_QuadsGrad0        = Kernel(program, "quads_grad0");
            m_Grad1             = Kernel(program, "grad1");

            m_Result     = cl::Buffer(m_Context, CL_MEM_WRITE_ONLY, sizeof(cl_scalar_t) );
            m_TempBuffer = cl::Buffer(m_Context, CL_MEM_READ_WRITE, sizeof(cl_scalar2_t) * max( 4*num_of_quads, 3*num_of_tris) );

            compute_delta();
        } catch ( Error &err )
        {
            std::cerr << "knupp_func :: error in " << err.what() << " code = " << err.err() << "\n";
        }
    }

    double rescale()
    {
        return 1.0;
    }

    Real eval_f(const Real &alpha = 0.0)
    {
        ++m_Stats.fevals;
        scalar_t result_quads = 0.0;
        scalar_t results_tris = 0.0;

        try{
            if( alpha != 0.0 )
            {
                m_LineSearchPrepare.setArg(0,m_num_of_inner);
                m_LineSearchPrepare.setArg(1,alpha);
                m_LineSearchPrepare.setArg(2,m_Points);
                m_LineSearchPrepare.setArg(3,m_Df);
                m_LineSearchPrepare.setArg(4,m_PointsAux);

                m_Queue.enqueueNDRangeKernel(
                    m_LineSearchPrepare,
                    cl::NullRange,
                    cl::NDRange(arch::clutils::round_up(m_num_of_threads, m_num_of_inner)),
                    cl::NDRange(m_num_of_threads));
            }

            if( m_num_of_triangles > 0 )
            {
                m_TrisFunc.setArg(0,m_num_of_triangles); 
                m_TrisFunc.setArg(1,m_D42); 
                m_TrisFunc.setArg(2,alpha != 0.0 ? m_PointsAux : m_Points);
                m_TrisFunc.setArg(3,m_Triangles);
                m_TrisFunc.setArg(4,m_TempBuffer);

                m_Queue.enqueueNDRangeKernel(
                    m_TrisFunc,
                    cl::NullRange,
                    cl::NDRange( arch::clutils::round_up(m_num_of_threads, m_num_of_triangles) ),
                    cl::NDRange(m_num_of_threads));
                
                arch::blas::avg( m_blas_hd, m_num_of_triangles, m_TempBuffer, &results_tris );
                results_tris = results_tris - 1.0;
            }

            if( m_num_of_quads > 0 ) 
            {
                m_QuadsFunc.setArg(0,m_num_of_quads); 
                m_QuadsFunc.setArg(1,m_D42); 
                m_QuadsFunc.setArg(2,alpha != 0.0 ? m_PointsAux : m_Points);
                m_QuadsFunc.setArg(3,m_Quads);
                m_QuadsFunc.setArg(4,m_TempBuffer);

                m_Queue.enqueueNDRangeKernel(
                    m_QuadsFunc,
                    cl::NullRange,
                    cl::NDRange( arch::clutils::round_up(m_num_of_threads,m_num_of_quads) ),
                    cl::NDRange(m_num_of_threads));
                  
                arch::blas::avg( m_blas_hd, m_num_of_quads, m_TempBuffer, &result_quads );
                result_quads = (result_quads/4.0) - 1.0;
            }
        } catch ( cl::Error &err ) {
            std::cerr << "eval_f :: error in " << err.what() << " code = " << err.err() << " error = " << arch::clutils::error(err.err() ) << "\n";
        }

        return results_tris + result_quads;
    }

    void eval_grad_f()
    {
        ++m_Stats.gradevals; 
        //compute_delta();

        try{
            if( m_num_of_triangles > 0 ) 
            {
                m_TrisGrad0.setArg(0,m_num_of_triangles); 
                m_TrisGrad0.setArg(1,m_D42);
                m_TrisGrad0.setArg(2,m_Points);
                m_TrisGrad0.setArg(3,m_Triangles);
                m_TrisGrad0.setArg(4,m_TempBuffer);

                m_Queue.enqueueNDRangeKernel(
                    m_TrisGrad0,
                    cl::NullRange,
                    cl::NDRange(round_up(m_num_of_threads,m_num_of_triangles) ),
                    cl::NDRange(m_num_of_threads));

                m_Grad1.setArg(0,m_num_of_inner ); 
                m_Grad1.setArg(1,m_Triangles_ia);
                m_Grad1.setArg(2,m_Triangles_ja);
                m_Grad1.setArg(3,m_TempBuffer);
                m_Grad1.setArg(4,1.0);
                m_Grad1.setArg(5,m_Df);

                m_Queue.enqueueNDRangeKernel(
                    m_Grad1,
                    cl::NullRange,
                    cl::NDRange(round_up(m_num_of_threads,m_num_of_inner) ),
                    cl::NDRange(m_num_of_threads));
            }

            if( m_num_of_quads > 0 ) 
            {
                m_QuadsGrad0.setArg(0,m_num_of_quads); 
                m_QuadsGrad0.setArg(1,m_D42);
                m_QuadsGrad0.setArg(2,m_Points);
                m_QuadsGrad0.setArg(3,m_Quads);
                m_QuadsGrad0.setArg(4,m_TempBuffer);

                m_Queue.enqueueNDRangeKernel(
                    m_QuadsGrad0,
                    cl::NullRange,
                    cl::NDRange(round_up(m_num_of_threads,m_num_of_quads) ),
                    cl::NDRange(m_num_of_threads));

                m_Grad1.setArg(0,m_num_of_inner ); 
                m_Grad1.setArg(1,m_Quads_ia);
                m_Grad1.setArg(2,m_Quads_ja);
                m_Grad1.setArg(3,m_TempBuffer);
                m_Grad1.setArg(4,1.0);
                m_Grad1.setArg(5,m_Df);

                m_Queue.enqueueNDRangeKernel(
                    m_Grad1,
                    cl::NullRange,
                    cl::NDRange(round_up(m_num_of_threads,m_num_of_inner) ),
                    cl::NDRange(m_num_of_threads));
            }

            //normalize the result
            cl_scalar_t scale_factor = 1.0 / ( 4.0*m_num_of_quads + m_num_of_triangles);
            arch::blas::scal( m_blas_hd, m_num_of_inner*2, scale_factor, m_Df, 1 );
        } catch ( cl::Error &err ) {
            std::cerr << "eval_grad_f :: error in " << err.what() << " code = " << err.err() << " error = " << arch::clutils::error(err.err() ) << "\n";
        }
    }

    void set_delta(const Real &delta)
    {
        m_D42 = 4.0 * delta * delta;
    }

    cl::Buffer &get_x() 
    {
        return m_Points;
    }

    cl::Buffer &get_df()
    {
        return m_Df;
    }

    struct stats
    {
        stats() : fevals(0), gradevals(0) {}
        cl_int fevals;
        cl_int gradevals;
    };

    stats get_stats()const { return m_Stats; }

    void check_derivatives(double eps = 1e-03)
    {
        std::cout << "Starting derivative checker...\n";
        std::vector< cl_scalar_t > px(2*m_num_of_inner);
        std::vector< cl_scalar_t > df(2*m_num_of_inner);
        std::vector< cl_scalar_t > x(2*m_num_of_inner);

        m_Queue.enqueueReadBuffer(m_Points, CL_TRUE, 0, x.size()*sizeof(cl_scalar_t), &x[0] );
        m_Queue.enqueueReadBuffer(m_Df, CL_TRUE, 0, df.size()*sizeof(cl_scalar_t), &df[0] );

        for( int i = 0; i < 2*m_num_of_inner; ++i ) {
            for( int j = 0; j < 2*m_num_of_inner; ++j )
                px[j] = ((i == j) ? (x[j]-1e-8) : x[j]);
            m_Queue.enqueueWriteBuffer(m_Points, CL_TRUE, 0, px.size()*sizeof(cl_scalar_t), &px[0] );
            cl_scalar_t f0 = eval_f();
            for( int j = 0; j < 2*m_num_of_inner; ++j )
                px[j] = ((i == j) ? (x[j]+1e-8) : x[j]);
            m_Queue.enqueueWriteBuffer(m_Points, CL_TRUE, 0, px.size()*sizeof(cl_scalar_t), &px[0] );
            cl_scalar_t f1 = eval_f();

            //cl_scalar_t scale_factor = 1.0 / ( 4.0*m_num_of_quads + m_num_of_triangles);
            cl_scalar_t fd = (f1 - f0) / (2*1e-08);

            if( fabs( fd - df[i] ) > eps ) {
                printf("grad[%d] = %.12lf~%.12lf [ %.12lf ]\n", i, fd, df[i], fabs( fd - df[i] ) );
                //std::cout << "grad_f[" << i << "] = " << fd << "~" << df[i] << "\n";
            }
        }
        std::cout << "derivative checker finished\n";

        m_Queue.enqueueWriteBuffer(m_Points, CL_TRUE, 0, x.size()*sizeof(cl_scalar_t), &x[0] );
    }

private:
    void compute_delta()
    {
        cl_scalar_t delta;

        try{
            int reduction_threads = m_num_of_threads;

            m_ComputeDelta.setArg(0,m_num_of_triangles); 
            m_ComputeDelta.setArg(1,m_num_of_quads);
            m_ComputeDelta.setArg(2,m_Points); 
            m_ComputeDelta.setArg(3,m_Triangles);
            m_ComputeDelta.setArg(4,m_Quads);
            m_ComputeDelta.setArg(5,sizeof(cl_scalar_t)*reduction_threads,NULL);
            m_ComputeDelta.setArg(6,m_Result);

            m_Queue.enqueueNDRangeKernel(
                m_ComputeDelta,
                cl::NullRange,
                cl::NDRange(reduction_threads),
                cl::NDRange(reduction_threads));

            m_Queue.enqueueReadBuffer(m_Result, CL_TRUE, 0, sizeof(scalar_t), &delta);
            if( delta != 0.0 ) std::cout << "D = " << delta << "\n";;
            set_delta(delta);

        } catch ( cl::Error &err ) {
            std::cerr << "compute_delta :: error in " << err.what() << " code = " << err.err() << " error = " << arch::clutils::error(err.err() ) << "\n";
        }
    }

    cl::Context m_Context;
    cl::Device m_Device;
    cl::CommandQueue &m_Queue;
    cl_int m_num_of_inner;
    cl_int m_num_of_triangles;
    cl_int m_num_of_quads;
    cl_int m_num_of_threads;
    cl_scalar_t m_D42;
    cl::Buffer &m_Points;
    cl::Buffer &m_PointsAux;
    cl::Buffer &m_Triangles;
    cl::Buffer &m_Triangles_ia;
    cl::Buffer &m_Triangles_ja;
    cl::Buffer &m_Quads;
    cl::Buffer &m_Quads_ia;
    cl::Buffer &m_Quads_ja;
    cl::Buffer &m_Df;
    cl::Buffer m_Result;
    cl::Buffer m_TempBuffer;
    cl::Kernel m_TrisFunc;
    cl::Kernel m_QuadsFunc;
    cl::Kernel m_ComputeDelta;
    cl::Kernel m_LineSearchPrepare;
    cl::Kernel m_TrisGrad0;
    cl::Kernel m_QuadsGrad0;
    cl::Kernel m_Grad1;
    arch::blas::handler m_blas_hd;
    stats m_Stats;
};

bool is_locked( vertex_ptr_t v )
{
    //Check if any edge of the vertex is free (boundary) or tjoin
    for( vertex_t::edge_iterator_t e = v->edges_begin();
        e != v->edges_end(); ++e )
    {
        if( is_free( *e ) || is_tjoin( *e ) )
            return true;
    }

    return !v->edges_size();        //There is the possibility that this is a free vertex
}

bool is_locked( face_ptr_t f )
{
    //Check if any edge of the vertex is free (boundary) or tjoin
    for( face_t::vertex_iterator_t v = f->verts_begin();
        v != f->verts_end(); ++v )
    {
        if( !is_locked(*v) && !(*v)->pinned )
            return false;
    }

    return true;
}

void build_csr( const cl::Context &context, const std::vector< int > &faces, int type, int offset, int inner_nodes, int max_neigbs, cl::Buffer &d_ia, cl::Buffer &d_ja)
{
    if( faces.empty() ) return;

#if 1
    std::vector< int > ia( inner_nodes );
    std::fill(ia.begin(),ia.end(),0);
    std::vector< int > ja( inner_nodes * max_neigbs );
    std::fill( ja.begin(), ja.end(), -1 );

    for( std::size_t i = 0; i < faces.size(); ++i ) {
        if( faces[i] >= 0 && faces[i] < inner_nodes ) {
            //std::cout << "j = " << (ia[faces[i]] * inner_nodes + faces[i]) << "\n";
            ja[ ia[faces[i]] * inner_nodes + faces[i] ] = int( (i / offset) * type + i%offset);
            ++ia[ faces[i] ];
        }
    }

    for( std::size_t i = 0; i < ia.size(); ++i ) 
        ia[i] = (int)(i + ia[i] * inner_nodes);
#else
    std::vector< int > ia(inner_nodes+1);
    std::fill(ia.begin(),ia.end(),0);
    
    for( std::size_t i = 0; i < faces.size(); ++i ){
        if( faces[i] >= 0 && faces[i] < inner_nodes ) {
            ++ia[faces[i]+1];
        }
    }
    
    for( std::size_t i = 0; i < inner_nodes; ++i )
        ia[i+1] += ia[i];
    
    std::vector< int > ja(ia.back());
    std::fill( ja.begin(), ja.end(), -1 );
    
    for( std::size_t i = 0; i < faces.size(); ++i ){
        if( faces[i] >= 0 && faces[i] < inner_nodes ) {
            for( std::size_t j = ia[faces[i]]; j < ia[faces[i]+1]; ++j )
                if( ja[j] == -1 ) {
                    ja[j] = (i / offset) * type + i%offset;
                    break;
                }
        }
    }
#endif

    d_ia = cl::Buffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(int) * ia.size(), &ia[0] );
    d_ja = cl::Buffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(int) * ja.size(), &ja[0] );


#if 0
    if( faces.empty() ) return;

    std::vector< int > ia(inner_nodes+1);
    std::fill(ia.begin(),ia.end(),0);

    for( std::size_t i = 0; i < faces.size(); ++i ){
        if( faces[i] >= 0 && faces[i] < inner_nodes ) {
            ++ia[faces[i]+1];
        }
    }

    for( std::size_t i = 0; i < inner_nodes; ++i )
        ia[i+1] += ia[i];

    std::vector< int > ja(ia.back());
    std::fill( ja.begin(), ja.end(), -1 );

    for( std::size_t i = 0; i < faces.size(); ++i ){
        if( faces[i] >= 0 && faces[i] < inner_nodes ) {
            for( std::size_t j = ia[faces[i]]; j < ia[faces[i]+1]; ++j )
                if( ja[j] == -1 ) {
                    ja[j] = (i / offset) * type + i%offset;
                    break;
                }
        }
    }

    d_ia = cl::Buffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(int) * ia.size(), &ia[0] );
    d_ja = cl::Buffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(int) * ja.size(), &ja[0] );
#endif
}

parameterization::SolverCL::SolverCL( 
    mesh_t &input_mesh,
    mesh_t &ouput_mesh,
    Options options) : m_Input(input_mesh), m_Output(ouput_mesh), m_Options(options)
{
    using namespace std;
    using namespace cl;
    using std::max;
    using std::min;

    try {
        int sel_platform = -1;
        int sel_device = -1;
        arch::clutils::init_best_device(sel_platform, sel_device, options.device_type);

        //arch::blas::profile(sel_platform, sel_device);

        std::vector<Platform> platforms;
        std::vector<cl::Device> devices;
        Platform::get(&platforms);

        cl_context_properties properties[] = 
        { CL_CONTEXT_PLATFORM, (cl_context_properties)(platforms[sel_platform])(), 0};

        m_Context = cl::Context(options.device_type, properties); 
        devices = m_Context.getInfo<CL_CONTEXT_DEVICES>();

        //create the command queue
        m_Queue = CommandQueue(m_Context, devices[sel_device], options.profile ? CL_QUEUE_PROFILING_ENABLE : 0);

        //Clean model : remove faces that are locked
        if( !options.free_boundaries ) {
            foreach( mesh_t::face_ptr_t f, input_mesh.faces() ) {
                if( is_locked(f) ) m_BannedFaces.push_back(f);
            }
        }

        std::cout << "Faces to erase = " << m_BannedFaces.size() << "\n";
        for( std::size_t i = 0; i < m_BannedFaces.size(); ++i ) 
            input_mesh.remove_face(m_BannedFaces[i],false);

        //first sort the vertices in this order : [free, fixed]
        m_verts_size = input_mesh.verts_size();
        m_fixed_pnts = 0;

        if( options.projection_type != 1 ) {
            for( std::size_t i = 0; i < m_verts_size - m_fixed_pnts; ) {
                mesh_t::vertex_ptr_t v = input_mesh.verts()[i];

                if( v->pinned || (!options.free_boundaries && is_locked( v )) ) {
                    input_mesh.swap_vertex( v, input_mesh.verts()[ m_verts_size - m_fixed_pnts - 1 ] );
                    ++m_fixed_pnts;
                } else {
                    ++i;
                }
            }
        }

        mesh_t::point_t pmin = input_mesh.verts()[0]->point();
        mesh_t::point_t pmax = input_mesh.verts()[0]->point();

        foreach( mesh_t::vertex_ptr_t v, input_mesh.verts() ) {
            mesh_t::point_t p = v->point();
            for( int i = 0; i < 3; i++ ) {
                pmin[i] = min(pmin[i],p[i]);
                pmax[i] = max(pmax[i],p[i]);
            }
        }

        m_Scale = input_mesh.edges()[0]->length();
        foreach( mesh_t::edge_ptr_t e, input_mesh.edges() ) {
            if( e->length() < m_Scale )
                m_Scale = e->length();
        }

        m_Scale *= 1e+02;
        m_Scale = 10.0;

        //m_Scale = arch::math::distance(pmax, pmin) * 0.01;
        m_Origin = pmin;
        m_Diagonal = arch::math::distance(pmax, pmin);
        std::size_t reverse_faces = 0;
        double min_area = 1000000.0;
      
        std::vector< cl_scalar_t > angles;
        std::vector< cl_scalar_t > factors;

        if( options.energy == 0 ) compute_mips_weights(angles, factors);

        if( options.projection_type == 1 ) {
            circular_projection(200000);

            if( !options.free_boundaries ) {
                for( std::size_t i = 0; i < m_verts_size - m_fixed_pnts; ) {
                    mesh_t::vertex_ptr_t v = input_mesh.verts()[i];

                    if( v->pinned  ) {
                        input_mesh.swap_vertex( v, input_mesh.verts()[ m_verts_size - m_fixed_pnts - 1 ] );
                        ++m_fixed_pnts;
                    } else {
                        ++i;
                    }
                }
            }

            foreach( mesh_t::vertex_ptr_t v, input_mesh.verts() ) {
                m_HostPnts.push_back( v->u * 10.0 );
                m_HostPnts.push_back( v->v * 10.0 );
            }
        } else if( options.projection_type == 0 ) {   //default mode (planar projection)
            //store the points in a buffer
            foreach( mesh_t::vertex_ptr_t v, input_mesh.verts() ) {
                if( !v->pinned ) {
                    m_HostPnts.push_back( (v->point().x - m_Origin.x) / m_Scale );
                    m_HostPnts.push_back( (v->point().y - m_Origin.y) / m_Scale );
                    //m_HostPnts.push_back( v->point().x / m_Scale );
                    //m_HostPnts.push_back( v->point().y / m_Scale );
                } else {
                    //m_HostPnts.push_back( (v->u - m_Origin.x) / m_Scale );
                    //m_HostPnts.push_back( (v->v - m_Origin.y) / m_Scale );
                    m_HostPnts.push_back( v->u );
                    m_HostPnts.push_back( v->v );
                }
            }
        } else {        //uv mode
            foreach( mesh_t::face_ptr_t f, input_mesh.faces() ) {
                std::vector< vec3d > proj;
                std::vector< vec3d > verts;
                foreach( face_t::vertex_ptr_t fv, f->verts() ) {
                    proj.push_back( vec3d( fv->u, fv->v, 0.0 ) );
                    verts.push_back( fv->point() );
                }

                vec3d e0 = proj[1] - proj[0];
                vec3d e1 = proj[2] - proj[1];
                vec3d e2 = proj[2] - proj[0];

                double det = e0.x * e2.y - e0.y * e2.x;
                if( det < 0.0 ) {
#if 0
                    std::cout << "Det = " << det << " face = " << f->id() << "\n";

                    std::cout << "v " << verts[0].x << " " << verts[0].y << " " << verts[0].z << "\n";
                    std::cout << "v " << verts[1].x << " " << verts[1].y << " " << verts[1].z << "\n";
                    std::cout << "v " << verts[2].x << " " << verts[2].y << " " << verts[2].z << "\n";
                  
                    std::cout << "vt " << proj[0].x << " " << proj[0].y << "\n";
                    std::cout << "vt " << proj[1].x << " " << proj[1].y << "\n";
                    std::cout << "vt " << proj[2].x << " " << proj[2].y << "\n";
#endif
                    reverse_faces++;
                    //break;
                }

                if( det < min_area ) 
                    min_area = det;

                //scalar_t det_w = max( e0.s0 * e2.s1 - e0.s1 * e2.s0, 1e-12 );
            }

            foreach( mesh_t::vertex_ptr_t v, input_mesh.verts() ) {
                m_HostPnts.push_back( v->u );
                m_HostPnts.push_back( v->v );
            }

#if 0
            //Energy contour plot
            double step = 0.025;
            foreach( mesh_t::vertex_ptr_t v, input_mesh.verts() ) {
                if( v->verts_size() > 3 ) {
                    double u_start = v->u - 0.75;
                    double v_start = v->v - 0.5;
                    std::vector< std::vector< double > > plot_2d;

                    for( double iu = u_start; iu < u_start + 1.25; iu += step ) {
                        plot_2d.push_back( std::vector< double >() );
                        for( double iv = v_start; iv < v_start + 1.25; iv += step ) {
                            v->u = iu;  v->v = iv;
                            double energy = 0.0;
                            foreach( mesh_t::face_ptr_t f, v->faces() ) {
                                std::vector< vec3d > proj;
                                foreach( face_t::vertex_ptr_t fv, f->verts() ) {
                                    proj.push_back( vec3d( fv->u, fv->v, 0.0 ) );
                                }

                                vec3d e0 = proj[1] - proj[0];
                                vec3d e1 = proj[2] - proj[1];
                                vec3d e2 = proj[2] - proj[0];
                                 
                                double det = e0.x * e2.y - e0.y * e2.x;
                                energy += det * det / powf( area( f->points() ), 2.0 ); 
                                //energy += det * det; 
                                
                            }
                            plot_2d.back().push_back( energy );
                            //std::cout << v->u << " " << v->v << " " << energy << "\n";
                        }
                    }

                    double iu = u_start;
                    for( std::size_t i = 0; i < plot_2d.size() - 1; ++i ) {
                        double iv = v_start;
                        for( std::size_t j = 0; j < plot_2d[i].size() - 1; ++j ) {
                             std::cout << iu << " " << iv << " " <<               plot_2d[i][j] << "\n";
                             std::cout << (iu+step) << " " << iv << " " <<        plot_2d[i+1][j] << "\n";
                             std::cout << (iu+step) << " " << (iv+step) << " " << plot_2d[i+1][j+1] << "\n";
                             std::cout << iu << " " << (iv+step) << " " <<        plot_2d[i][j+1] << "\n\n";
                    
                            iv += step;
                        }

                        std::cout << "\n\n";
                       
                        iu += step;
                    }
                }
            }
#endif
        }

        bool reverse_order = (2*reverse_faces) >= input_mesh.faces_size();
        std::cout << "Model scaling = " << m_Scale << "\n";
        std::cout << "Reverse faces = " << reverse_faces << "/" << input_mesh.faces_size() << "\n";
        std::cout << "Reverse order = "  << reverse_order << "\n";
        std::cout << "Min area = " << min_area << "\n";
        std::cout << "Model verts = " << m_verts_size << "\n";
        std::cout << "Model fixed verts = " << m_fixed_pnts << "\n";

        if( !angles.empty() )
            m_TrisAngles = cl::Buffer( m_Context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(cl_double)*angles.size(), &angles[0] );

        if( !factors.empty() )
            m_TrisFactors = cl::Buffer( m_Context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(cl_double)*factors.size(), &factors[0] );

        //store the faces in a buffer
        std::vector< int > tris;
        std::vector< int > quads;

        foreach( mesh_t::face_ptr_t f, input_mesh.faces() ) {
            if( f->verts_size() == 4 ) {
                foreach( face_t::vertex_ptr_t fv, f->verts() ) {
                    quads.push_back( int(fv->id()) );
                }
            } else if( f->verts_size() == 3 ) {
                foreach( face_t::vertex_ptr_t fv, f->verts() ) {
                    tris.push_back( int(fv->id()) );
                }
                if( reverse_order ) std::swap( tris.back(), tris[ tris.size() - 3] );
                tris.push_back( -1 );
            }
        }
        int max_valence = 0;
        foreach( mesh_t::vertex_ptr_t v, input_mesh.verts() ) 
            max_valence = int(max( v->verts_size(), max_valence ));

        m_tris_size  = tris.size()  / 4;
        m_quads_size = quads.size() / 4;

        std::cout << "Triangles = " << m_tris_size << "\n";
        std::cout << "Quads = " << m_quads_size << "\n";
        std::cout << "Max valence = " << max_valence << "\n";
        std::size_t inner_verts = m_verts_size - m_fixed_pnts;

        //build ia,ja matrics
        build_csr( m_Context, tris, 3, 4, inner_verts, max_valence, m_Tris_ia, m_Tris_ja );
        build_csr( m_Context, quads, 4, 4, inner_verts, max_valence, m_Quads_ia, m_Quads_ja );
       
        m_Pnts     = cl::Buffer( m_Context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(cl_scalar_t)*m_HostPnts.size(), &m_HostPnts[0] );
        m_PntsAux  = cl::Buffer( m_Context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(cl_scalar_t)*m_HostPnts.size(), &m_HostPnts[0] );

        if( m_tris_size > 0 ) 
            m_Tris  = cl::Buffer( m_Context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(cl_int)*tris.size()     , &tris[0]);
        else 
            m_Tris  = cl::Buffer( m_Context, CL_MEM_READ_ONLY, sizeof(cl_int) );

        if( m_quads_size > 0 )
            m_Quads = cl::Buffer( m_Context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(cl_int)*quads.size()    , &quads[0]);
        else
            m_Quads = cl::Buffer( m_Context, CL_MEM_READ_ONLY , sizeof(cl_int) );

        m_Df    = cl::Buffer( m_Context, CL_MEM_READ_WRITE, sizeof(cl_scalar2_t)*inner_verts );

    } catch ( Error &err ) {
        std::cerr << "SolverCL :: Error in : " << err.what() << " code = " << err.err() << " error = " << arch::clutils::error(err.err() ) << "\n";
    }
}

parameterization::SolverCL::~SolverCL()
{

}

bool parameterization::SolverCL::solve(parameterization::Stats &solve_stats)
{
    using std::min;

    //solve_stats.elapsed_ms = timeGetTime();
    solve_stats.success = true;
    solve_stats.iterations = 0;

    if( m_Options.energy == 0 ) 
    {
        mips_func<cl_scalar_t> kf(
            m_Queue, 
            m_verts_size,
            m_verts_size - m_fixed_pnts, 
            m_tris_size, 
            m_quads_size, 
            m_Pnts, 
            m_PntsAux, 
            m_Tris, 
            m_Tris_ia, 
            m_Tris_ja, 
            m_TrisAngles,
            m_TrisFactors,
            m_Quads,
            m_Quads_ia,
            m_Quads_ja,
            m_Df,
            m_Options.theta);
     
        if( 1 )
        {
#if 0
            std::vector< float > HostPntsFloat( m_HostPnts.size() );
            for( int i = 0; i < m_HostPnts.size(); ++i )
                HostPntsFloat[i] = m_HostPnts[i];

            cl::Buffer new_buffer = cl::Buffer( m_Context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(cl_float)*HostPntsFloat.size(), &HostPntsFloat[0] );
            std::swap( m_Pnts, new_buffer );
#endif

            solve_stats.elapsed_ms = timeGetTime();
            //Laplace iterations
            double mix_factor = 0.0;
            double f = 0.0, f_old = 0.0;
            std::cout.precision(16);
#if 0
            std::cout << "Smoothing\n";
            double smooth_factor = 0.5;
            f_old = kf.eval_laplace_f(smooth_factor, true );
            for( int i = 1; i < m_Options.un_iters+1; ++i ) 
            {
                //std::cout << i << "\n";
                f = kf.eval_laplace_f(smooth_factor, (i % 1000 == 0) );
                
                if( i % 1000 == 0 ) {
                    std::cout << i << " " << f << " " << smooth_factor << " " << (f_old - f) << "\n";
                    //if( f_old - f < 1e-04 ) {
                    if( f_old - f < 1e-12 ) {
                        //if( smooth_factor > 0.1 && f_old != f  ) 
                        if( smooth_factor > 1e-12 && f_old != f  ) 
                            smooth_factor *= 0.5;
                        else
                            break;
                    }

                    f_old = f;
                }
            }
#endif
#if 1

            std::cout << "Untangling\n";
            float un_factor = 0.5f;
            f_old = f = kf.eval_untangle_f(un_factor, true );
            for( int i = 1; i < m_Options.un_iters+1; ++i ) 
            {
                //std::cout << i << "\n";
                f = kf.eval_untangle_f(un_factor, (i % 1000 == 0) );
                //f = kf.eval_untangle_f(un_factor, true );
                //f = kf.eval_untangle_f(un_factor, true );
                //std::cout << i << " " << f << " " << un_factor << " " << (f_old - f) << "\n";
               
                if( i % 1000 == 0 ) {
                    std::cout << i << " " << f << " " << un_factor << " " << (f_old - f) << "\n";
                    if( f_old - f <= 1e-12 ) {
                        if( un_factor > 0.00001 && f_old != f ) {
                            un_factor *= 0.5;
                        } else break;
                    }
                
                    f_old = f;
                }
            }
#endif
#if 0
            m_Queue.enqueueReadBuffer(m_Pnts, CL_TRUE, 0, HostPntsFloat.size()*sizeof(cl_float), &HostPntsFloat[0] );
            for( int i = 0; i < HostPntsFloat.size(); ++i ) 
                m_HostPnts[i] = HostPntsFloat[i];

            std::swap( m_Pnts, new_buffer );
            m_Queue.enqueueWriteBuffer(m_Pnts   , CL_TRUE, 0, m_HostPnts.size() * sizeof(cl_scalar_t), &m_HostPnts[0]);
            m_Queue.enqueueWriteBuffer(m_PntsAux, CL_TRUE, 0, m_HostPnts.size() * sizeof(cl_scalar_t), &m_HostPnts[0]);
#endif

            //Recompute delta to see if the mesh is untangled
            kf.compute_delta();
        }

        if( m_Options.opt_iters ) {
            if( m_Options.scale_iters < (~(0u) >> 1) )
                kf.rescale();  

            //kf.eval_grad_f();
            //kf.check_derivatives();
            //std::cout << "f = " << kf.eval_f() << "\n";
            line_search_func_t func = boost::bind( &mips_func<cl_scalar_t>::eval_f, &kf, _1 );
            base_line_search<cl_scalar_t, line_search_func_t> kls( func );
            solve_stats.elapsed_ms = timeGetTime();
            
            //if( m_Options.scale_iters < (~(0u) >> 1) )
            //    kf.rescale();  

            //kf.eval_grad_f();
            //for( double a = -1.5; a <= 5.0; a += 0.001 ) 
            //    std::cout << a << " " << kf.eval_f(a) << "\n";
            //
            //arch::blas::handler hd(m_Queue);
            //arch::blas::print( hd,2, kf.get_df() ); 

            //arch::blas::axpy( hd, 2, -0.11, kf.get_df(), 1, kf.get_x(), 1);
            //std::cout << "Final = " << kf.eval_f() << "\n";
            //kf.rescale();       //Initial problem scaling
            
            
            arch::solver::pcg<1,cl_scalar_t>(m_Queue, kf, kls, 2*(m_verts_size - m_fixed_pnts), m_Options.opt_iters, m_Options.scale_iters ); 
        }

        m_Scale = kf.get_total_scale();
    } else {
        knupp_func<cl_scalar_t> kf(
            m_Queue, 
            m_verts_size - m_fixed_pnts, 
            m_tris_size, 
            m_quads_size, 
            m_Pnts, 
            m_PntsAux, 
            m_Tris, 
            m_Tris_ia, 
            m_Tris_ja, 
            m_Quads,
            m_Quads_ia,
            m_Quads_ja,
            m_Df);
        
        line_search_func_t func = boost::bind( &knupp_func<cl_scalar_t>::eval_f, &kf, _1 );
        base_line_search<cl_scalar_t, line_search_func_t> kls( func );
        solve_stats.elapsed_ms = timeGetTime();
       
        arch::solver::pcg<1,cl_scalar_t>(m_Queue, kf, kls, 2*(m_verts_size - m_fixed_pnts), m_Options.opt_iters, m_Options.scale_iters ); 
    }
   
    m_Queue.enqueueReadBuffer(m_Pnts, CL_TRUE, 0, m_HostPnts.size()*sizeof(cl_scalar_t), &m_HostPnts[0] );
   
    solve_stats.elapsed_ms = timeGetTime() - solve_stats.elapsed_ms;

    //Reinsert the banned faces (Faces with all the vertices locked that were removed during the optimization process)
    for( std::size_t i = 0; i < m_BannedFaces.size(); ++i ) 
        m_Input.add_face(m_BannedFaces[i]);

    for( std::size_t i = 0; i < m_HostPnts.size(); i += 2 )
    {
        //std::cout << i << "/" << m_HostPnts.size() << " while " << m_Input.verts_size() << "\n";
        //uid_t id = m_Output.add_vertex( vertex_ptr_t( new vertex_t(
        //    m_HostPnts[i],m_HostPnts[i+1],0.0 ) ) );
        //uid_t id = m_Output.add_vertex( vertex_ptr_t( new vertex_t(
        //    m_HostPnts[i]*m_Scale + m_Origin.x,m_HostPnts[i+1]*m_Scale + m_Origin.y,0.0 ) ) );
        uid_t id = m_Output.add_vertex( vertex_ptr_t( new vertex_t(
            m_Input.verts()[i/2]->point() ) ) );
        
       //m_Output.verts()[id]->u = (m_HostPnts[i]*m_Scale + m_Origin.x) / m_Diagonal;
       //m_Output.verts()[id]->v = (m_HostPnts[i+1]*m_Scale + m_Origin.y) / m_Diagonal;
       //m_Output.verts()[id]->u = m_HostPnts[i] / 10.0;
       //m_Output.verts()[id]->v = m_HostPnts[i+1] / 10.0; 
       
        m_Output.verts()[id]->u = m_HostPnts[i]   / m_Scale; 
        m_Output.verts()[id]->v = m_HostPnts[i+1] / m_Scale; 

        m_Output.verts()[id]->pinned = m_Input.verts()[i/2]->pinned;
    }

    for( mesh_t::face_iterator_t f = m_Input.faces_begin(); f != m_Input.faces_end(); ++f )
    {
        std::vector< vertex_ptr_t > verts;
        std::vector< vertex_ptr_t > indices((*f)->verts_begin(),(*f)->verts_end());

        for( std::size_t i = 0 ; i < indices.size(); ++i )
            verts.push_back( m_Output.verts()[ indices[i]->id() ] );

        m_Output.add_face( face_ptr_t( new face_t(verts.begin(),verts.end()) ) );
    }
    
    return true;
}

void parameterization::SolverCL::compute_mips_weights(std::vector< cl_scalar_t > &angles, std::vector< cl_scalar_t > &factors)
{
    //double mean_value = 0.0;
    //std::cout << "compute_mips_weights()\n";
    foreach( mesh_t::face_ptr_t f, m_Input.faces() )
    {
        if( !f->is_triangle() ) continue;

        std::vector< point_t > pnts(f->points_begin(),f->points_end());

        cl_scalar_t a0 = acos( dot( 
            normalize( pnts[1] - pnts[0] ), 
            normalize( pnts[2] - pnts[0] ) ) );

        cl_scalar_t a1 = acos( dot( 
            normalize( pnts[2] - pnts[1] ), 
            normalize( pnts[0] - pnts[1] ) ) );

        cl_scalar_t a2 = acos( dot( 
            normalize( pnts[0] - pnts[2] ), 
            normalize( pnts[1] - pnts[2] ) ) );

         //Clamp the weights between 5 and 85 degrees
         //if( a0 < 0.087266462599716474f ) a0 = 0.087266462599716474f;
         //else if( a0 > 1.4835298641951802f ) a0 = 1.4835298641951802f;
         //
         //if( a1 < 0.087266462599716474f ) a1 = 0.087266462599716474f;
         //else if( a1 > 1.4835298641951802f ) a1 = 1.4835298641951802f;
         //
         //if( a2 < 0.087266462599716474f ) a2 = 0.087266462599716474f;
         //else if( a2 > 1.4835298641951802f ) a2 = 1.4835298641951802f;

        angles.push_back( 1.0 / tan(a0) );
        angles.push_back( 1.0 / tan(a1) );
        angles.push_back( 1.0 / tan(a2) );
        //angles.push_back( arch::math::distance_sq( pnts[1], pnts[0]) );
        angles.push_back( 2.0*area( f->points() ) );

        scalar_t cotA = 1.0 / tan(a0);
        scalar_t cotB = 1.0 / tan(a1);
        
        //std::vector< scalar_t > ang;
        //ang.push_back( 1.0/ tan(a0) );
        //ang.push_back( 1.0/ tan(a1) );
        //ang.push_back( 1.0/ tan(a2) );
        //std::sort( ang.begin(), ang.end() );

        
#if 0
        //printf("%f : %f\n", cotA + cotB, ang[1] + ang[2] );
        //mean_value += ang[1] + ang[2];
        factors.push_back( cotA );
        factors.push_back( cotB );
        factors.push_back( 1.0 / tan(a2) );
        factors.push_back( 2.0*area( f->points() ) );
        //factors.push_back( 1.0 / (2.0 * (cotA + cotB) * area( f->points() )) );
        //factors.push_back( cotA + cotB );
#endif

        double scale = 1.0 / sqrt( 2.0 * area( f->points() ) * (cotA + cotB) );
        factors.push_back( scale );
        factors.push_back( scale * cotA );
        factors.push_back( scale * cotB );

        //factors.push_back( 0.0 );
    }

}

bool find_boundary( mesh_t::edge_ptr_t e, std::vector< mesh_t::vertex_ptr_t > &boundary, double &perimeter ) 
{
    boundary.clear();
    perimeter = 0.0;
    std::map< mesh_t::edge_ptr_t, bool > checked;
    std::deque< mesh_t::edge_ptr_t > q;
    mesh_t::edge_ptr_t cur_edge = e;
    mesh_t::vertex_ptr_t cur_vert = e->v0();
    mesh_t::vertex_ptr_t prev_vert = cur_vert;

    while( 1 ) {
        bool found_edge = false;
        boundary.push_back( cur_vert );
        
        perimeter += cur_edge->length();
        
        //mark the edge as checked
        checked[cur_edge] = true;
        prev_vert = cur_vert;
        cur_vert = prev_vert == cur_edge->v0() ? cur_edge->v1() : cur_edge->v0();
        foreach( mesh_t::edge_ptr_t ee, cur_vert->edges() ) {
            if( is_free( ee ) && ee != cur_edge ) {
                //found closed boundary
                if( ee == e ) return true;

                //do not select checked edge
                if( checked.find(ee) != checked.end() ) continue;

                found_edge = true;
                cur_edge = ee;
                break;
            }
        }

        //failed to find closed boundary
        if( !found_edge ) 
            return false;
    }
}

void parameterization::SolverCL::untangle_laplace(int max_iters)
{
#if 0
    using std::min;
    using std::max;
    using namespace cl;

    std::size_t verts_size = m_Input.verts_size();
    std::size_t fixed_pnts = 0;

    for( std::size_t i = 0; i < m_verts_size - m_fixed_pnts; ) {
        mesh_t::vertex_ptr_t v = m_Input.verts()[i];

        if( v->pinned || (!m_Options.free_boundaries && is_locked( v )) ) {
            m_Input.swap_vertex( v, m_Input.verts()[ m_verts_size - m_fixed_pnts - 1 ] );
            ++m_fixed_pnts;
        } else {
            ++i;
        }
    }
    
    //Read the source file
    ifstream source_file("laplace.cl");
    std::string source_code(
        std::istreambuf_iterator<char>(source_file),
        (std::istreambuf_iterator<char>()));

    Program::Sources source(1, std::make_pair(source_code.c_str(), source_code.size()+1));
    Program program = Program(m_Context, source);

    cl::Device device =  m_Queue.getInfo<CL_QUEUE_DEVICE>();	
    std::vector< Device > devices;
    devices.push_back(device);

    char options[255];
    try{
        sprintf(options, "-DUSE_DOUBLE=%d  -cl-mad-enable -cl-fast-relaxed-math",CL_USE_DOUBLE);
        program.build(devices,options);
    } catch ( Error &err )
    {
        std::cerr << err.what() << "\n";
        std::cerr << program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(device) << "\n";
    }

    cl::Kernel laplace_kernel = Kernel(program, "laplace_kernel");

    //compute the vertex data
    std::vector< cl_scalar_t > x0;
    
    foreach( mesh_t::vertex_ptr_t v, m_Input.verts() ) {
        x0.push_back( v->u );
        x0.push_back( v->v );
    }

    std::vector< int > ia(1,0);
    std::vector< int > ja;
    //compute the neighboring data
    foreach( mesh_t::vertex_ptr_t v, m_Input.verts() ) {
        ia.push_back( ia.back() + v->verts_size() );
        foreach( mesh_t::vertex_ptr_t vv, v->ordered_verts() ) {
        }

        //foreach( mesh_t::vertex_ptr_t vv, v->verts() ) {
        //    ja.push_back( vv->id() );
        //}
    }

    //copy data to the device
    try{
        cl::Buffer d_x0( m_Context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(cl_scalar_t) * x0.size(), &x0[0] );
        cl::Buffer d_x1( m_Context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(cl_scalar_t) * x0.size(), &x0[0] );
        cl::Buffer d_ia( m_Context, CL_MEM_READ_ONLY  | CL_MEM_COPY_HOST_PTR, sizeof(cl_int)      * ia.size(), &ia[0] );
        cl::Buffer d_ja( m_Context, CL_MEM_READ_ONLY  | CL_MEM_COPY_HOST_PTR, sizeof(cl_int)      * ja.size(), &ja[0] );
        cl_int inner_verts = verts_size - fixed_pnts;

        laplace_kernel.setArg(0,inner_verts); 
        laplace_kernel.setArg(2,d_ia);
        laplace_kernel.setArg(3,d_ja);

        //execute the kernel
        for( int k = 0; k < max_iters; ++k ) {
            if( k % 1000 == 0 ) 
                std::cout << k << "/" << max_iters << "\n";

            laplace_kernel.setArg(1,d_x0);
            laplace_kernel.setArg(4,d_x1);

            m_Queue.enqueueNDRangeKernel(
                laplace_kernel,
                cl::NullRange,
                cl::NDRange(inner_verts),
                cl::NullRange);

            std::swap( d_x0, d_x1 );
        }

        //read back the data
        m_Queue.enqueueReadBuffer(d_x0, CL_TRUE, 0, x0.size()*sizeof(cl_scalar_t), &x0[0] );
    } catch ( cl::Error &err ) {
        std::cerr << "circular_projection :: error in " << err.what() << " code = " << err.err() << " error = " << arch::clutils::error(err.err() ) << "\n";
    }

    //copy data back to the mesh
    foreach( mesh_t::vertex_ptr_t v, m_Input.verts() ) {
        v->u = (x0[v->id()*2 + 0] / radius);
        v->v = (x0[v->id()*2 + 1] / radius);
    }
#endif
}


void parameterization::SolverCL::circular_projection(int max_iters)
{
    using std::min;
    using std::max;
    using namespace cl;

    std::vector< mesh_t::vertex_ptr_t > boundary;
    double perimeter = 0;

    //find the boundary with the longest perimeter
    foreach( mesh_t::edge_ptr_t e, m_Input.edges() ) {
        double p;
        std::vector< mesh_t::vertex_ptr_t > b;

        if( is_free( e ) && find_boundary( e, b, p ) ) {
            if( p > perimeter ) {
                perimeter = p;
                boundary = b;
            }
        }
    }

    mesh_t::point_t pmin = m_Input.verts()[0]->point();
    mesh_t::point_t pmax = m_Input.verts()[0]->point();
    mesh_t::point_t c(0.0);

    foreach( mesh_t::vertex_ptr_t v, m_Input.verts() ) {
        mesh_t::point_t p = v->point();
        for( int i = 0; i < 3; i++ ) {
            pmin[i] = min(pmin[i],p[i]);
            pmax[i] = max(pmax[i],p[i]);
        }

        c += v->point();
    }

    c /= m_Input.verts_size();

    double angle = 0.0;
    double PI_2 = 8.0 * atan(1.0);
    foreach( mesh_t::vertex_ptr_t v, m_Input.verts() ) {
        v->u = v->point().x - c.x;
        v->v = v->point().y - c.y;
        //v->u = v->v = 0.0;
    }

    double radius = max( pmax[0] - pmin[0], pmax[1] - pmin[1] );
    std::cout << "Radius = " << radius << "\n";

    for( std::size_t i = 0; i < boundary.size(); ++i ) {
        angle += PI_2 * (arch::math::distance( boundary[(i+1)%boundary.size()]->point() , boundary[i]->point() ) / perimeter );
        boundary[i]->u = cos( -angle ) * radius;
        boundary[i]->v = sin( -angle ) * radius;
        boundary[i]->pinned = true;
    }

    std::size_t verts_size = m_Input.verts_size();
    std::size_t fixed_pnts = 0;

    //put the boundary vertices at the end
    for( std::size_t i = 0; i < verts_size - fixed_pnts; ) {
        mesh_t::vertex_ptr_t v = m_Input.verts()[i];

        if( v->pinned ) {
            m_Input.swap_vertex( v, m_Input.verts()[ verts_size - fixed_pnts - 1 ] );
            ++fixed_pnts;
        } else {
            ++i;
        }
    }
    
    //Read the source file
    ifstream source_file("laplace.cl");
    std::string source_code(
        std::istreambuf_iterator<char>(source_file),
        (std::istreambuf_iterator<char>()));

    Program::Sources source(1, std::make_pair(source_code.c_str(), source_code.size()+1));
    Program program = Program(m_Context, source);

    cl::Device device =  m_Queue.getInfo<CL_QUEUE_DEVICE>();	
    std::vector< Device > devices;
    devices.push_back(device);

    char options[255];
    try{
        sprintf(options, "-DUSE_DOUBLE=%d  -cl-mad-enable -cl-fast-relaxed-math",CL_USE_DOUBLE);
        program.build(devices,options);
    } catch ( Error &err )
    {
        std::cerr << err.what() << "\n";
        std::cerr << program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(device) << "\n";
    }

    cl::Kernel laplace_kernel = Kernel(program, "laplace_kernel");

    //compute the vertex data
    std::vector< cl_scalar_t > x0;
    
    foreach( mesh_t::vertex_ptr_t v, m_Input.verts() ) {
        x0.push_back( v->u );
        x0.push_back( v->v );
    }

    std::vector< int > ia(1,0);
    std::vector< int > ja;
    //compute the neighboring data
    foreach( mesh_t::vertex_ptr_t v, m_Input.verts() ) {
        ia.push_back( ia.back() + v->verts_size() );
        foreach( mesh_t::vertex_ptr_t vv, v->verts() ) {
            ja.push_back( vv->id() );
        }
    }

    //copy data to the device
    try{
        cl::Buffer d_x0( m_Context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(cl_scalar_t) * x0.size(), &x0[0] );
        cl::Buffer d_x1( m_Context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(cl_scalar_t) * x0.size(), &x0[0] );
        cl::Buffer d_ia( m_Context, CL_MEM_READ_ONLY  | CL_MEM_COPY_HOST_PTR, sizeof(cl_int)      * ia.size(), &ia[0] );
        cl::Buffer d_ja( m_Context, CL_MEM_READ_ONLY  | CL_MEM_COPY_HOST_PTR, sizeof(cl_int)      * ja.size(), &ja[0] );
        cl_int inner_verts = verts_size - fixed_pnts;

        laplace_kernel.setArg(0,inner_verts); 
        laplace_kernel.setArg(2,d_ia);
        laplace_kernel.setArg(3,d_ja);

        //execute the kernel
        for( int k = 0; k < max_iters; ++k ) {
            if( k % 1000 == 0 ) 
                std::cout << k << "/" << max_iters << "\n";

            laplace_kernel.setArg(1,d_x0);
            laplace_kernel.setArg(4,d_x1);

            m_Queue.enqueueNDRangeKernel(
                laplace_kernel,
                cl::NullRange,
                cl::NDRange(inner_verts),
                cl::NullRange);

            std::swap( d_x0, d_x1 );
        }

        //read back the data
        m_Queue.enqueueReadBuffer(d_x0, CL_TRUE, 0, x0.size()*sizeof(cl_scalar_t), &x0[0] );
    } catch ( cl::Error &err ) {
        std::cerr << "circular_projection :: error in " << err.what() << " code = " << err.err() << " error = " << arch::clutils::error(err.err() ) << "\n";
    }

    //copy data back to the mesh
    foreach( mesh_t::vertex_ptr_t v, m_Input.verts() ) {
        v->u = (x0[v->id()*2 + 0] / radius);
        v->v = (x0[v->id()*2 + 1] / radius);
    }
}

void parameterization::SolverCL::compute_weights()
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

            w = (1.0f / tan(a0) + 1.0f / tan(a1));
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
