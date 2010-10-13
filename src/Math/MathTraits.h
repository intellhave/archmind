/*
  Archmind Non-manifold Geometric Kernel
  Copyright (C) 2010 Athanasiadis Theodoros

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

#ifndef MATH_TRAITS_H
#define MATH_TRAITS_H

#include <cmath>
#include <cfloat>

namespace math
{

template<typename Real>
struct traits
{
	static const Real pi;
	static const Real zero_tol;
	static const Real epsilon;
	static const Real max_real;
};

template<typename Real> const Real traits<Real>::pi(4.0f*atan(4.0f));
template<typename Real> const Real traits<Real>::zero_tol(1e-06f);
template<typename Real> const Real traits<Real>::epsilon(FLT_EPSILON);
template<typename Real> const Real traits<Real>::max_real(FLT_MAX);

//inline template<> const double traits<double>::pi(4.0*atan(4.0));
//inline template<> const double traits<double>::zero_tol(1e-08);
//inline template<> const double traits<double>::epsilon(DBL_EPSILON);
//inline template<> const double traits<double>::max_real(DBL_MAX);

};

#endif
