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
class vector2
{
public:
    //typedef Real value_type;
    typedef Real real_t;

    // constructors
    vector2<Real>( const vector2<Real> &v );
    vector2<Real>( const Real &x_, const Real &y_ );
    vector2<Real>( const Real &v );
    vector2<Real>( void );

    // operators
    void operator+=( const vector2<Real> &v );
    void operator-=( const vector2<Real> &v);
    void operator/=(const Real &v);
    void operator*=(const Real &v);

    vector2<Real> operator-();

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
Real dot(const vector2<Real> &a, const vector2<Real> &b);

//!Absolute dot Product of two vectors
template< typename Real>
Real adot(const vector2<Real> &a, const vector2<Real> &b);

template< typename Real>
Real magnitude(const vector2<Real> &a);

template<typename Real>
std::ostream &operator<<(std::ostream &output, const vector2<Real> &v);

template< typename Real >
class vector3
{
public:
    typedef Real real_t;

    template< typename T >
    vector3<Real>( const vector3<T> &v ) : 
        x(v.x), y(v.y), z(v.z)
    { 
        x = Real(v.x); y = Real(v.y); z = Real(v.z);
    }

    template< typename T >
    vector3<Real>( const T &x_, const Real &y_, const Real &z_) :
        x(x_), y(y_), z(z_)
    {
    }
   
    vector3<Real>( void );
    vector3<Real>( const Real *ptr ); 
    vector3<Real>( const Real &x_, const Real &y_, const Real &z_);
    vector3<Real>( const Real &v );

    // operators
    bool operator==(const vector3<Real> &v);

    vector3<Real> &operator=( const Real *ptr ); 

    void operator+=( const vector3<Real> &v );
    void operator-=( const vector3<Real> &v);
    void operator/=(const Real &v);
    void operator*=(const Real &v);
    
    vector3<Real> operator-()const;

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
Real dot(const vector3<Real> &a, const vector3<Real> &b);

//!Absolute dot Product of two vectors
template< typename Real>
Real adot(const vector3<Real> &a, const vector3<Real> &b);

//!magnitude of vector
template< typename Real>
Real magnitude(const vector3<Real> &a);

//!cross Product of two vectors
template< typename Real>
vector3<Real> cross(const vector3<Real> &a, const vector3<Real> &b);

//!normal Vector of three co-planar vectors (e.g. a triangle)
template< typename Real>
vector3<Real> normal(const vector3<Real> &v1, const vector3<Real> &v2, const vector3<Real> &v3);

template<typename Real>
std::ostream &operator<<(std::ostream &output, const vector3<Real> &v);

template< typename Real >
class vector4
{
public:
    //typedef Real value_type;
    typedef Real real_t;
    
    // constructors
    vector4<Real>( const vector4<Real> &v );
    vector4<Real>( const Real &x_, const Real &y_, const Real &z_, const Real &w_);
    vector4<Real>( const Real &v);
    vector4<Real>( void );

    // operators
    void operator+=( const vector4<Real> &v );
    void operator-=( const vector4<Real> &v );
    void operator/=(const Real &v );
    void operator*=(const Real &v );

    vector4<Real> operator-();
    
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

//!normalize vector
template< typename VectorType >
VectorType normalize(const VectorType &a);

template< typename Real >
Real clip_line_to_plane
(const vector3<Real> &lp, const vector3<Real> &ld, const vector3<Real> &pn, const Real &pd)
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
(const vector3<Real> &point, const vector3<Real> &v1, const vector3<Real> &v2, const vector3<Real> &v3)
{
    vector3<Real> E0, E1, Q;
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

typedef vector2<double> vector2d;
typedef vector2<float> vector2f;
typedef vector2<int> vector2i;
typedef vector2<unsigned int> vector2u;
typedef vector2<short> vector2s;
typedef vector2<unsigned short> vector2us;

typedef vector3<double> vector3d;
typedef vector3<float> vector3f;
typedef vector3<int> vector3i;
typedef vector3<unsigned int> vector3u;
typedef vector3<short> vector3s;
typedef vector3<unsigned short> vector3us;
typedef vector3<unsigned char> vector3ub;

typedef vector4<double> vector4d;
typedef vector4<float> vector4f;
typedef vector4<int> vector4i;
typedef vector4<unsigned int> vector4u;
typedef vector4<short> vector4s;
typedef vector4<unsigned short> vector4us;
typedef vector4<unsigned char> vector4ub;

#include "Vector.inl"

}

}

#endif
