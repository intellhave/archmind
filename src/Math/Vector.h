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

#ifndef MATH_VECTOR_H
#define MATH_VECTOR_H

#include "MathTraits.h"
#include <cstdlib>      //for abs
#include <cmath>        //for sqrt
#include <iostream>

namespace arch
{

namespace math
{

template< typename Real >
class vec2
{
public:
    //typedef Real value_type;
    typedef Real real_t;

    // constructors
    vec2<Real>( const vec2<Real> &v );
    vec2<Real>( const Real &x_, const Real &y_ );
    vec2<Real>( const Real &v );
    vec2<Real>( void );

    // operators
    void operator+=( const vec2<Real> &v );
    void operator-=( const vec2<Real> &v);
    void operator/=(const Real &v);
    void operator*=(const Real &v);

    vec2<Real> operator-();

    Real &operator()(unsigned int index); 
    Real operator()(unsigned int index)const;

    Real &operator[](unsigned int index);
    Real operator[](unsigned int index)const;
    
    std::size_t size()const;

    // variables
    Real x, y;
};

//!dot Product of two vectors
template< typename Real>
Real dot(const vec2<Real> &a, const vec2<Real> &b);

//!Absolute dot Product of two vectors
template< typename Real>
Real adot(const vec2<Real> &a, const vec2<Real> &b);

template< typename Real>
Real magnitude(const vec2<Real> &a);

template<typename Real>
std::ostream &operator<<(std::ostream &output, const vec2<Real> &v);

template< typename Real >
class vec3
{
public:
    typedef Real real_t;

    template< typename T >
    vec3<Real>( const vec3<T> &v ) : 
        x(v.x), y(v.y), z(v.z)
    { 
        x = Real(v.x); y = Real(v.y); z = Real(v.z);
    }

    template< typename T >
    vec3<Real>( const T &x_, const Real &y_, const Real &z_) :
        x(x_), y(y_), z(z_)
    {
    }
   
    vec3<Real>( void );
    vec3<Real>( const Real *ptr ); 
    vec3<Real>( const Real &x_, const Real &y_, const Real &z_);
    vec3<Real>( const Real &v );

    // operators
    bool operator==(const vec3<Real> &v);

    vec3<Real> &operator=( const Real *ptr ); 

    void operator+=( const vec3<Real> &v );
    void operator-=( const vec3<Real> &v);
    void operator/=(const Real &v);
    void operator*=(const Real &v);
    
    vec3<Real> operator-()const;

    Real &operator()(unsigned int index);
    Real operator()(unsigned int index)const;

    Real &operator[](unsigned int index);
    Real operator[](unsigned int index)const;
    
    std::size_t size()const;

    // variables
    Real x, y, z;
};

//!dot Product of two vectors
template< typename Real>
Real dot(const vec3<Real> &a, const vec3<Real> &b);

//!Absolute dot Product of two vectors
template< typename Real>
Real adot(const vec3<Real> &a, const vec3<Real> &b);

//!magnitude of vector
template< typename Real>
Real magnitude(const vec3<Real> &a);

//!cross Product of two vectors
template< typename Real>
vec3<Real> cross(const vec3<Real> &a, const vec3<Real> &b);

//!normal Vector of three co-planar vectors (e.g. a triangle)
template< typename Real>
vec3<Real> normal(const vec3<Real> &v1, const vec3<Real> &v2, const vec3<Real> &v3);

template<typename Real>
std::ostream &operator<<(std::ostream &output, const vec3<Real> &v);

template< typename Real >
class vec4
{
public:
    //typedef Real value_type;
    typedef Real real_t;
    
    // constructors
    vec4<Real>( const vec4<Real> &v );
    vec4<Real>( const Real &x_, const Real &y_, const Real &z_, const Real &w_);
    vec4<Real>( const Real &v);
    vec4<Real>( void );

    // operators
    void operator+=( const vec4<Real> &v );
    void operator-=( const vec4<Real> &v );
    void operator/=(const Real &v );
    void operator*=(const Real &v );

    vec4<Real> operator-();
    
    Real &operator()(unsigned int index);
    Real operator()(unsigned int index)const;

    Real &operator[](unsigned int index);
    Real operator[](unsigned int index)const;
    
    std::size_t size()const;

    // variables
    Real x, y, z, w;
};

//!distance of two vectors
template< typename VectorType >
typename VectorType::real_t distance(const VectorType &a, const VectorType &b);

//!squared distance of two vectors
template< typename VectorType >
typename VectorType::real_t distance_sq(const VectorType &a, const VectorType &b);

//!normalize vector
template< typename VectorType >
VectorType normalize(const VectorType &a);

template< typename Real >
Real clip_line_to_plane
(const vec3<Real> &lp, const vec3<Real> &ld, const vec3<Real> &pn, const Real &pd)
{
    using std::abs;

    Real dotND = dot(pn,ld);

    //check if the line is parallel to the plane
    if( abs( dotND ) < traits<Real>::zero_tol )
        return traits<Real>::max_real;

    return (+pd -dot(lp,pn)) / dotND;
}

template< typename Real >
bool point_in_triangle
(const vec3<Real> &point, const vec3<Real> &v1, const vec3<Real> &v2, const vec3<Real> &v3)
{
    vec3<Real> E0, E1, Q;
    Real e00, e01, e11, Delta, q0, q1, s0, s1;

    E0 = v2 - v1;
    E1 = v3 - v1;
    
    Q = point - v1;

    e00 = dot(E0,E0);
    e01 = dot(E0,E1);
    e11 = dot(E1,E1);

    q0 = dot(E0,Q);
    q1 = dot(E1,Q);

    s0 = e11 * q0 - e01 * q1;
    s1 = e00 * q1 - e01 * q0;
    Delta = e00 * e11 - (e01 * e01);

    if( (s0 >= -traits<Real>::zero_tol && s1 >= -traits<Real>::zero_tol) && ((s0 + s1) / Delta) <= (1.0 + traits<Real>::zero_tol) )
        return true;
    else
        return false;
}

typedef vec2<double> vec2d;
typedef vec2<float> vec2f;
typedef vec2<int> vec2i;
typedef vec2<unsigned int> vec2u;
typedef vec2<short> vec2s;
typedef vec2<unsigned short> vec2us;

typedef vec3<double> vec3d;
typedef vec3<float> vec3f;
typedef vec3<int> vec3i;
typedef vec3<unsigned int> vec3u;
typedef vec3<short> vec3s;
typedef vec3<unsigned short> vec3us;
typedef vec3<unsigned char> vec3ub;

typedef vec4<double> vec4d;
typedef vec4<float> vec4f;
typedef vec4<int> vec4i;
typedef vec4<unsigned int> vec4u;
typedef vec4<short> vec4s;
typedef vec4<unsigned short> vec4us;
typedef vec4<unsigned char> vec4ub;

#include "Vector.inl"

}

}

#endif
