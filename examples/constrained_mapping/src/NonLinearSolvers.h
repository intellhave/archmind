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

#ifndef OPENCL_NON_LINEAR_SOLVER_H
#define OPENCL_NON_LINEAR_SOLVER_H

#define __CL_ENABLE_EXCEPTIONS
#include <CL/cl.hpp>
#include "OpenCLBlas.h"
#include <iostream>


#ifndef CL_USE_DOUBLE
#define CL_USE_DOUBLE 1
#endif


namespace arch
{
  namespace solver
  {
    template<typename Real, typename Func,  typename LineSearch>
    bool cg(cl::CommandQueue &q, Func &f, LineSearch &ls, const int n, const int max_iters = 50)
    {
      using std::max;
      using std::min;

      //1:d(0) = r(0) = -f'(x0)
      //while not converged
      //2:find a(i) that minimizes f(x(i)+a(i)*d(i))
      //3:x(i+1)=x(i)+a(i)*d(i)
      //4:r(i+1)=-f'(x(i+1))
      //5:b=max{dot(r(i+1),r(i+1)-r(i)),0)
      //6:d(i+1)=r(i+1)+b*d(i)
      cl::Context context = q.getInfo<CL_QUEUE_CONTEXT>();
      cl::Device device = q.getInfo<CL_QUEUE_DEVICE>();

      cl::Buffer r0(context, CL_MEM_READ_WRITE, sizeof(Real)*n );
      cl::Buffer d1(context, CL_MEM_READ_WRITE, sizeof(Real)*n );
      cl::Buffer &x = f.get_x();
      cl::Buffer &r1 = f.get_df();

      blas::handler h(q);		

      //1:d(0) = r(0) = -f'x(0)
      f.eval_grad_f();
      //blas::print(h, n, r1,true);
      blas::copy(h, n, r1, 1, r0, 1); 

      Real alpha, beta = 0.0;
      Real dot_r0r0, dot_r1rn;
      Real fmin;
      Real wolfe_acc, wolfe_c1 = 0.5;

      alpha = -1.0;

      for( int iter = 1; iter < max_iters; ++iter )
      {
        //2:find a(i) that minimizes f(x(i)+a(i)*d(i))
        //fmin = ls.line_search(lower_limit,0.0,alpha,1e-6,256);
        fmin = ls.line_search(alpha,1e-6,8);
        
        //3:x(i+1)=x(i)+a(i)*d(i)
        blas::axpy(h, n, alpha, r1, 1, x, 1);

        //store d(i)
        blas::copy(h, n, r1, 1, d1, 1);

        //4:r(i+1)=-f'(x(i+1))
        f.eval_grad_f();

        //compute dot( r(i), r(i) )
        blas::dot(h, n, r0, 1, r0, 1, &dot_r0r0);

        //compute r(i) = -r(i)
        blas::scal(h, n, -1.0, r0, 1);

        //compute r(i) = r(i+1) - r(i)
        blas::axpy(h, n, 1.0, r1, 1, r0, 1);

        //compute dot( r(i+1), r(i+1) - r(i) )
        blas::dot(h, n, r1, 1, r0, 1, &dot_r1rn);

        //5:b=max{dot(r(i+1),r(i+1)-r(i)/dot(r(i),r(i)),0)
        beta = max(dot_r1rn/dot_r0r0, 0.0);
        if( dot_r0r0 < 1e-08 || beta > 1e+04 || (iter % n == 0) ) 
          beta = 0.0;

        if( iter == 1 || iter % 100 == 0 ) 
        std::cout << iter << " " << fmin << " " << f.get_stats().fevals << " " << alpha << " " << beta << "\n";

        //std::cout << "iter = " << iter << ", fmin = " << fmin << ", alpha = " << alpha << ", beta = " << beta << "\n";

        //copy r(i) = r(i+1)
        blas::copy(h, n, r1, 1, r0, 1);

        //r(i+1) = r(i+1) + b*d(i)
        blas::axpy(h, n, beta, d1, 1, r1, 1);
      }

      return true;
    }

    template<int m, typename Real, typename Func, typename LineSearch>
    bool pcg(cl::CommandQueue &queue, Func &f, LineSearch &ls, int n, int max_iters = 50, int scale_iters = ((unsigned)~0) >> 1)
    {
      using std::max;
      using std::min;
      using std::abs;

      std::cout << "Preconditioned cg with m = " << m << "\n";

      cl::Context context = queue.getInfo<CL_QUEUE_CONTEXT>();
      cl::Device device = queue.getInfo<CL_QUEUE_DEVICE>();

      Real a[m] = {0.0};
      Real dot_ys[m] = {0.0};
      Real hscale[m] = {1.0};
      cl::Buffer s[m];
      for( int i = 0; i < m; ++i )
        s[i] = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(Real)*n );
      cl::Buffer y[m];
      for( int i = 0; i < m; ++i )
        y[i] = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(Real)*n );

      cl::Buffer &x = f.get_x();
      cl::Buffer xp(context, CL_MEM_READ_WRITE, sizeof(Real)*n );       //holds the previous x
      cl::Buffer &d = f.get_df();                                       //holds the search direction and the current gradient 
      cl::Buffer gp(context, CL_MEM_READ_WRITE, sizeof(Real)*n );       //holds the previous gradient
      cl::Buffer dp(context, CL_MEM_READ_WRITE, sizeof(Real)*n );       //holds the previous search direction
      
      blas::handler h(queue);		
      //1:d(0) = r(0) = -f'x(0)
      f.eval_grad_f();

      //Real min_der, max_der;
      //blas::amin(h, n, d, &min_der);
      //blas::amax(h, n, d, &max_der);
      //
      //std::cout << "Min derivative = " << min_der << "\n";
      //std::cout << "Max derivative = " << max_der << "\n";
      //f.check_derivatives();
      
      Real alpha, b, beta;
      Real fmin = 0.0;
      Real pfmin;

      //a0 = a1 = a2 = alpha = -1.0;
      Real dot_yy, dot_yg, dot_yd;
      //Real wolfe_acc, wolfe_c1 = 1e-01;
      Real wolfe_acc, wolfe_c1 = 0.3;

      blas::copy(h, n, d, 1, gp, 1); 

      //std::cout << f.eval_f() << "\n";
      
      int iter = 1;
      int end = 0;
      int k = 1;
     
      for( ; k < max_iters; ++k )
      {
        //copy previous vectors
        blas::copy(h, n, d, 1, dp, 1);
        blas::copy(h, n, x, 1, xp, 1);

        blas::dot(h, n, gp, 1, d, 1, &wolfe_acc);
        wolfe_acc *= wolfe_c1;

        //find a(i) that minimizes f(x(i)+a(i)*d(i))
        pfmin = fmin;

        //The initial decent direction may have a poor scaling therefore we perform the initial line search with higher accuracy
        if( iter == 1 ) fmin = ls.line_search(alpha,1e-16,30, wolfe_acc);
        else fmin = ls.line_search(alpha,1e-6,8, wolfe_acc);
        
        //fmin = ls.line_search(alpha,1e-16,30, 10.0);
      
        //fmin = ls.line_search(alpha,1e-15,256);
        //fmin = ls.line_search(alpha,1e-12,30, wolfe_acc);
        //if( (abs( fmin - ofmin ) / max( abs(fmin), abs(fmin) )) < 1e-06 ) break;
        //std::cout << (abs( fmin - ofmin ) / max( abs(fmin), abs(fmin) )) << "\n";
        //fmin = ls.line_search(lower_limit,0.0,alpha,1e-13,256);

        if( k % scale_iters == 0 ) {
            double scale_factor = f.rescale();
            iter = 1; end = 0;
            blas::copy(h, n, gp, 1, d, 1);
            ls.restart();
            if( fabs( scale_factor - 1.0 ) < 1e-03 )    //Disable scaling
                scale_iters = max_iters;

            continue;
        }

        if( alpha >= -1e-16 ) 
        {
          if( iter == 1 ) { 
              std::cout << "Failed to find a decent direction (alpha = " << alpha << ")\n";
              break;
          }

          iter = 1; end = 0;
          blas::copy(h, n, gp, 1, d, 1);
          ls.restart();
          std::cout << "cg restart\n";
          continue;
        }

        //if( k % 25 == 0 && pfmin - fmin < 1e-08 ) {
        //    std::cout << "Covergence (" << pfmin - fmin << ")\n";
        //    break;
        //}

        //x(i+1)=x(i)+a(i)*d(i)
        blas::axpy(h, n, alpha, d, 1, x, 1);
        
        f.eval_grad_f();
       
        //blas::copy(h, n, x, 1, s[end], 1);           //store x_{k+1}
        //blas::axpy(h, n, -1.0, xp, 1, s[end], 1);    //s_k = x_{k+1} - x_k
        blas::xmy(h, n, x, xp, s[end]);         //s_k = x_{k+1} - x_k
        //blas::copy(h, n, d, 1, y[end], 1);           //store g_k
        //blas::axpy(h, n, -1.0, gp, 1, y[end], 1);    //y_k = g_{k+1} - g_k
        blas::xmy(h, n, d, gp, y[end]);         //y_k = g_{k+1} - g_k
        
        blas::copy(h, n, d, 1, gp, 1); 
       
        blas::dot(h, n, y[end], 1, d, 1, &dot_yg );
        blas::dot(h, n, y[end], 1, dp, 1, &dot_yd );
       
        beta = dot_yg/dot_yd;
        //std::cout << dot_yd << "\n";
        if( dot_yd < 1e-12 || beta > 1e+04 || (iter % n == 0) ) 
          beta = 0.0;
      
        blas::dot(h, n, y[end], 1, s[end], 1, dot_ys+end);
        blas::dot(h, n, y[end], 1, y[end], 1, &dot_yy);

        hscale[end] = dot_ys[end] / dot_yy;
        Real scale = max( hscale[(end+1)%m] , hscale[end] );
       
        int bound = (m <= iter) ? m : iter;
        ++iter;

        end = (end + 1) % m;
        int j = end;
        for( int i = 0; i < bound; ++i ) 
        {
          j = (j + m - 1) % m;  //cyclic shift
          blas::dot( h, n, s[j], 1, d, 1, a+j );
          a[j] /= dot_ys[j];
          blas::axpy(h, n, -a[j], y[j], 1, d, 1 );
        }

        //z = H^0q, here we can apply scaling to the basic hessian
        blas::scal(h, n, scale, d, 1);

        for( int i = 0; i < bound; ++i ) 
        {
          blas::dot(h, n, y[j], 1, d, 1, &b);
          b /= dot_ys[j];
          blas::axpy(h, n, a[j] - b, s[j], 1, d, 1);
          j = (j + 1) % m;
        }

        blas::axpy(h, n, beta, dp, 1, d, 1 );

        if( k == 1 || k % 100 == 0 ) 
        std::cout << k << " " << fmin << " " << f.get_stats().fevals << " " << alpha << " " << beta << " " << scale << "\n";
      }
      std::cout << k << " " << fmin << " " << f.get_stats().fevals << " " << alpha << " " << beta << " " << "\n";

      return true;
    }

   template<typename Real, typename Func, typename LineSearch>
    bool hpcg(cl::CommandQueue &queue, Func &f, LineSearch &ls, const int num_of_m, int n, const int max_iters = 50)
    {
      using std::max;
      using std::min;
      using std::abs;
      int m = num_of_m;

      std::cout << "Preconditioned cg with m = " << m << "\n";

      cl::Context context = queue.getInfo<CL_QUEUE_CONTEXT>();
      cl::Device device = queue.getInfo<CL_QUEUE_DEVICE>();

      std::vector< Real > a(m, 0.0);
      std::vector< Real > dot_ys(m, 0.0);
      std::vector< Real > hscale(m, 1.0);

      std::vector< cl::Buffer > s(m);
      for( int i = 0; i < m; ++i )
        s[i] = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(Real)*n );
      std::vector< cl::Buffer > y(m);
      for( int i = 0; i < m; ++i )
        y[i] = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(Real)*n );

      cl::Buffer &x = f.get_x();
      cl::Buffer xp(context, CL_MEM_READ_WRITE, sizeof(Real)*n );       //holds the previous x
      cl::Buffer &d = f.get_df();                                       //holds the search direction and the current gradient 
      cl::Buffer gp(context, CL_MEM_READ_WRITE, sizeof(Real)*n );       //holds the previous gradient
      cl::Buffer dp(context, CL_MEM_READ_WRITE, sizeof(Real)*n );       //holds the previous search direction
      
      blas::handler h(queue);		
      //1:d(0) = r(0) = -f'x(0)
      f.eval_grad_f();
 
      //f.check_derivatives();
      
      Real alpha, b, beta;
      Real fmin;
      Real pfmin;

      //a0 = a1 = a2 = alpha = -1.0;
      Real dot_yy, dot_yg, dot_yd;
      //Real wolfe_acc, wolfe_c1 = 1e-01;
      Real wolfe_acc, wolfe_c1 = 0.5;

      blas::copy(h, n, d, 1, gp, 1); 

      //std::cout << f.eval_f() << "\n";
      
      int iter = 1;
      int end = 0;
      int k = 1;
      for( ; k < max_iters; ++k )
      {
        //copy previous vectors
        blas::copy(h, n, d, 1, dp, 1);
        blas::copy(h, n, x, 1, xp, 1);

        blas::dot(h, n, gp, 1, d, 1, &wolfe_acc);
        wolfe_acc *= wolfe_c1;

        //find a(i) that minimizes f(x(i)+a(i)*d(i))
        pfmin = fmin;
        //std::cout << "Line search\n";
        fmin = ls.line_search(alpha,1e-06,8, wolfe_acc);
        //fmin = ls.line_search(alpha,1e-15,256);
        //fmin = ls.line_search(alpha,1e-12,30, wolfe_acc);
        //if( (abs( fmin - ofmin ) / max( abs(fmin), abs(fmin) )) < 1e-06 ) break;
        //std::cout << (abs( fmin - ofmin ) / max( abs(fmin), abs(fmin) )) << "\n";
        //fmin = ls.line_search(lower_limit,0.0,alpha,1e-13,256);
        
        if( k > 1 && (k % 1000 == 0) && m > 1 ) { 
            m = m / 2;
            std::cout << "m is now " << m << "\n";
            iter = 1; end = 0;
            blas::copy(h, n, gp, 1, d, 1);
            ls.restart();
            //std::cout << "cg restart\n";
            continue;
        }

        if( alpha >= 0.0 ) 
        {
          if( iter == 1 ) { 
              std::cout << "Failed to find a decent direction\n";
              break;
          }

          iter = 1; end = 0;
          blas::copy(h, n, gp, 1, d, 1);
          ls.restart();
          std::cout << "cg restart\n";
          continue;
        }

        if( k % 25 == 0 && pfmin - fmin < 1e-08 ) break;

        //x(i+1)=x(i)+a(i)*d(i)
        blas::axpy(h, n, alpha, d, 1, x, 1);
        
        f.eval_grad_f();
       
        blas::copy(h, n, x, 1, s[end], 1);           //store x_{k+1}
        blas::axpy(h, n, -1.0, xp, 1, s[end], 1);    //s_k = x_{k+1} - x_k

        blas::copy(h, n, d, 1, y[end], 1);        //store g_k
        blas::axpy(h, n, -1.0, gp, 1, y[end], 1); //y_k = g_{k+1} - g_k

        blas::copy(h, n, d, 1, gp, 1); 
       
        blas::dot(h, n, y[end], 1, d, 1, &dot_yg );
        blas::dot(h, n, y[end], 1, dp, 1, &dot_yd );
       
        beta = dot_yg/dot_yd;
        //std::cout << dot_yd << "\n";
        if( dot_yd < 1e-12 || beta > 1e+04 || (iter % n == 0) ) 
          beta = 0.0;
      
        blas::dot(h, n, y[end], 1, s[end], 1, &dot_ys[end]);
        blas::dot(h, n, y[end], 1, y[end], 1, &dot_yy);

        hscale[end] = dot_ys[end] / dot_yy;
        Real scale = max( hscale[(end+1)%m] , hscale[end] );
       
        int bound = (m <= iter) ? m : iter;
        ++iter;

        end = (end + 1) % m;
        int j = end;
        for( int i = 0; i < bound; ++i ) 
        {
          j = (j + m - 1) % m;  //cyclic shift
          blas::dot( h, n, s[j], 1, d, 1, &a[j] );
          a[j] /= dot_ys[j];
          blas::axpy(h, n, -a[j], y[j], 1, d, 1 );
        }

        //z = H^0q, here we can apply scaling to the basic hessian
        blas::scal(h, n, scale, d, 1);

        for( int i = 0; i < bound; ++i ) 
        {
          blas::dot(h, n, y[j], 1, d, 1, &b);
          b /= dot_ys[j];
          blas::axpy(h, n, a[j] - b, s[j], 1, d, 1);
          j = (j + 1) % m;
        }

        blas::axpy(h, n, beta, dp, 1, d, 1 );

        if( k == 1 || k % 100 == 0 ) 
        std::cout << k << " " << fmin << " " << f.get_stats().fevals << " " << alpha << " " << beta << " " << scale << "\n";
      }
      std::cout << k << " " << fmin << " " << f.get_stats().fevals << " " << alpha << " " << beta << " " << "\n";

      return true;
    }
  }

};


#endif
