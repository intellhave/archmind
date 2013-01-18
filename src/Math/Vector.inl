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

/*
 * vec2<Real> Class Constructors 
 */
template< typename Real >    
inline vec2<Real>::vec2( const vec2<Real> &v )
{
    x = v.x; y = v.y;
}

template< typename Real >    
inline vec2<Real>::vec2( const Real &x_, const Real &y_ )
{
    x = x_; y = y_;
}

template< typename Real >    
inline vec2<Real>::vec2( const Real &v ) 
{
    x = v; y = v;
}

template< typename Real >    
inline vec2<Real>::vec2( void )
{}
      
/*
 * vec2<Real> Member Operators
 */

template< typename Real >    
inline void vec2<Real>::operator+=( const vec2<Real> &v )
{
    x += v.x;
    y += v.y;
}

template< typename Real >    
inline void vec2<Real>::operator-=( const vec2<Real> &v )
{
    x -= v.x;
    y -= v.y;
}

template< typename Real >    
inline void vec2<Real>::operator/=(const Real &v)
{
    x /= v;
    y /= v;
}

template< typename Real >    
inline void vec2<Real>::operator*=(const Real &v)
{
    x *= v;
    y *= v;
}

template< typename Real >    
inline vec2<Real> vec2<Real>::operator-()
{
    const vec2<Real> v(-x,-y);
    return v;
}

template< typename Real >    
inline Real &vec2<Real>::operator()(unsigned int index)
{
    return *(((Real *) &x) + index); 
}

template< typename Real >    
inline Real vec2<Real>::operator()(unsigned int index)const
{ 
    return *(((Real *) &x) + index); 
}

template< typename Real >    
inline Real &vec2<Real>::operator[](unsigned int index)
{
    return *(((Real *) &x) + index); 
}

template< typename Real >    
inline Real vec2<Real>::operator[](unsigned int index)const
{ 
    return *(((Real *) &x) + index); 
}

template< typename Real >
std::size_t vec2<Real>::size()const
{
    return 2;
}

/*
 * vec2<Real> Non Member Operators
 */
template< typename Real >
inline vec2<Real> operator+( const vec2<Real> &a, const vec2<Real> &b )
{
    const vec2<Real> v( a.x+b.x,a.y+b.y );
    return v;
}

template< typename Real >
inline vec2<Real> operator-( const vec2<Real> &a, const vec2<Real> &b )
{
    const vec2<Real> v(a.x-b.x,a.y-b.y);
    return v;
}

template< typename Real >
inline vec2<Real> operator/( const vec2<Real> &a, const Real b)
{
    const vec2<Real> v(a.x/b,a.y/b);
    return v;
}

template< typename Real >
inline vec2<Real> operator*( const vec2<Real> &a, const Real b)
{
    const vec2<Real> v(a.x*b,a.y*b);
    return v;
}

template< typename Real >
inline vec2<Real> operator*( const Real &b, const vec2<Real> &a)
{
    const vec2<Real> v(a.x*b,a.y*b);
    return v;
}

template< typename Real >
inline bool operator==(const vec2<Real> &a, const vec2<Real> &b )
{
    return (a.x == b.y) && (a.y == b.y);
}

/* 
 * vec2<Real> Functions
 */
template< typename Real>
inline Real magnitude(const vec2<Real> &a)
{
    return sqrt(a.x*a.x + a.y*a.y);
}

template< typename Real>
inline Real dot(const vec2<Real> &a, const vec2<Real> &b)
{
    return (a.x * b.x + a.y * b.y);
}

template< typename Real>
inline Real adot(const vec2<Real> &a, const vec2<Real> &b)
{
    using std::abs;
    return abs(dot(a,b));
}

template<typename Real>
inline std::ostream &operator << (std::ostream &output, const vec2<Real> &v)
{
    output << v.x << ',' << v.y;
    
    return output;
}


/*
 * vec3<Real> Class Constructors
 */


template< typename Real >
inline vec3<Real>::vec3( const Real *ptr )
{
    x = ptr[0]; y = ptr[1]; z = ptr[2];
}

template< typename Real >
inline vec3<Real>::vec3( const Real &x_, const Real &y_, const Real &z_)
{
    x = x_; y = y_; z = z_;
}

template< typename Real >
inline vec3<Real>::vec3( const Real &v ) 
{
    x = v; y = v; z = v;
}

template< typename Real >
inline vec3<Real>::vec3( void )
{}

/*
 * vec3<Real> Member Operators
 */

template< typename Real >
inline bool vec3<Real>::operator==(const vec3<Real> &v )
{
    using std::abs;

#if 0
    return
        std::abs(x - v.x) < (Real)0.0001 &&
        std::abs(y - v.y) < (Real)0.0001 &&
        std::abs(z - v.z) < (Real)0.0001;
#endif

    return x == v.x && y == v.y && z == v.z;
}

template< typename Real >
inline void vec3<Real>::operator+=( const vec3<Real> &v )
{
    x += v.x;
    y += v.y;
    z += v.z;
}

template< typename Real >
inline void vec3<Real>::operator-=( const vec3<Real> &v)
{
    x -= v.x;
    y -= v.y;
    z -= v.z;
}

template< typename Real >
inline void vec3<Real>::operator/=(const Real &v)
{
    x /= v;
    y /= v;
    z /= v;
}

template< typename Real >
inline void vec3<Real>::operator*=(const Real &v)
{
    x *= v;
    y *= v;
    z *= v;
}

template< typename Real >
inline vec3<Real> vec3<Real>::operator-()const
{
    const vec3<Real> v(-x,-y,-z);
    return v;
}

template< typename Real >
inline Real &vec3<Real>::operator()(unsigned int index)
{ 
    return *(((Real *) &x) + index); 
}

template< typename Real >
inline Real vec3<Real>::operator()(unsigned int index)const
{
    return *(((Real *) &x) + index); 
}

template< typename Real >
inline Real &vec3<Real>::operator[](unsigned int index)
{
    return *(((Real *) &x) + index); 
}

template< typename Real >
inline Real vec3<Real>::operator[](unsigned int index)const
{ 
    return *(((Real *) &x) + index); 
}

template< typename Real >
std::size_t vec3<Real>::size()const
{
    return 3;
}

/*
 * vec3<Real> Non Member Operators
 */
template< typename Real >
inline vec3<Real> operator+( const vec3<Real> &a, const vec3<Real> &b )
{
       const vec3<Real> v( a.x+b.x,a.y+b.y, a.z+b.z );
       return v;
}

template< typename Real >
inline vec3<Real> operator-( const vec3<Real> &a, const vec3<Real> &b )
{
       const vec3<Real> v(a.x-b.x,a.y-b.y,a.z-b.z);
       return v;
}


template< typename Real >
inline vec3<Real> operator/( const vec3<Real> &a, const Real &f )
{
       const vec3<Real> v(a.x/f,a.y/f,a.z/f);
       return v;
}

template< typename Real >
inline vec3<Real> operator*( const vec3<Real> &a, const Real &f )
{
       const vec3<Real> v(a.x*f,a.y*f,a.z*f);
       return v;
}

template< typename Real >
inline vec3<Real> operator*( const Real &f, const vec3<Real> &a )
{
       const vec3<Real> v(a.x*f,a.y*f,a.z*f);
       return v;
}

template< typename Real >
inline bool operator==(const vec3<Real> &a, const vec3<Real> &b )
{
       return (a.x == b.y) && (a.y == b.y) && (a.z == b.z);
}

template<typename Real>
inline std::ostream &operator << (std::ostream &output, const vec3<Real> &v)
{
    //output << v.x << ',' << v.y << ',' << v.z;
    output << v.x << ' ' << v.y << ' ' << v.z;
    
    return output;
}

/*
 * vec3<Real> Functions
 */
template< typename Real>
inline Real dot(const vec3<Real> &a, const vec3<Real> &b)
{
    return (a.x * b.x + a.y * b.y + a.z * b.z);
}

template< typename Real>
inline Real adot(const vec3<Real> &a, const vec3<Real> &b)
{
    using std::abs;
    return abs(dot(a,b));
}


template< typename Real>
inline Real magnitude(const vec3<Real> &a)
{
    using std::sqrt;

    return sqrt(a.x*a.x + a.y*a.y + a.z*a.z);
}

template< typename Real>
inline vec3<Real> cross(const vec3<Real> &a, const vec3<Real> &b)
{
    const vec3<Real> vcross(
            ((a.y * b.z) - (a.z * b.y)),
            ((a.z * b.x) - (a.x * b.z)),
            ((a.x * b.y) - (a.y * b.x)));

    return vcross;
}

template< typename Real>
inline vec3<Real> normal(const vec3<Real> &v1, const vec3<Real> &v2, const vec3<Real> &v3)
{
    return normalize(cross(v3-v1,v2-v1));   
}

/*
 * vec4<Real> Class Constructors 
 */
template< typename Real >
inline vec4<Real>::vec4( const vec4<Real> &v )
{
    x = v.x; y = v.y; z = v.z; w = v.w;
}
          
template< typename Real >
inline vec4<Real>::vec4( const Real &x_, const Real &y_, const Real &z_, const Real &w_)
{
    x = x_; y = y_; z = z_; w = w_;
}
          
template< typename Real >
inline vec4<Real>::vec4( const Real &v ) 
{
    x = v; y = v; z = v; w = v;
}
          
template< typename Real >
inline vec4<Real>::vec4( void )
{}

/*
 * vec4<Real> Member Operators
 */

template< typename Real >
inline void vec4<Real>::operator+=( const vec4<Real> &v )
{
    x += v.x; y += v.y; z += v.z; w += v.w;
}

template< typename Real >
inline void vec4<Real>::operator-=( const vec4<Real> &v)
{
    x -= v.x; y -= v.y; z -= v.z; w -= v.w;
}

template< typename Real >
inline void vec4<Real>::operator/=(const Real &v)
{
    x /= v; y /= v; z /= v; w /= v;
}

template< typename Real >
inline void vec4<Real>::operator*=(const Real &v)
{
    x *= v; y *= v; z *= v; z *= v;
}

template< typename Real >
inline vec4<Real> vec4<Real>::operator-()
{
    const vec4<Real> v(-x,-y,-z,-w);
    return v;
}


template< typename Real >
inline Real &vec4<Real>::operator()(unsigned int index)
{ 
    return *(((Real *) &x) + index); 
}
          
template< typename Real >
inline Real vec4<Real>::operator()(unsigned int index)const
{ 
    return *(((Real *) &x) + index); 
}

template< typename Real >
inline Real &vec4<Real>::operator[](unsigned int index)
{ 
    return *(((Real *) &x) + index); 
}
          
template< typename Real >
inline Real vec4<Real>::operator[](unsigned int index)const
{ 
    return *(((Real *) &x) + index); 
}

template< typename Real >
std::size_t vec4<Real>::size()const
{
    return 4;
}

/*
 * vec4<Real> Non Member Operators
 */
template< typename Real >
inline vec4<Real> operator+( const vec4<Real> &a, const vec4<Real> &b )
{
       const vec4<Real> v( a.x+b.x,a.y+b.y, a.z+b.z, a.w+b.w );
       return v;
}

template< typename Real >
inline vec4<Real> operator-( const vec4<Real> &a, const vec4<Real> &b )
{
       const vec4<Real> v(a.x-b.x,a.y-b.y,a.z-b.z,a.w-b.w);
       return v;
}

template< typename Real >
inline vec4<Real> operator/( const vec4<Real> &a, const Real &f )
{
       const vec4<Real> v(a.x/f,a.y/f,a.z/f,a.w/f);
       return v;
}

template< typename Real >
inline vec4<Real> operator*( const vec4<Real> &a, const Real &f )
{
       const vec4<Real> v(a.x*f,a.y*f,a.z*f,a.w*f);
       return v;
}

template< typename Real >
inline vec4<Real> operator*( const Real f, const vec4<Real> &a )
{
       const vec4<Real> v(a.x*f,a.y*f,a.z*f,a.w*f);
       return v;
}

template< typename Real >
inline bool operator==(const vec4<Real> &a, const vec4<Real> &b )
{
       return (a.x == b.y) && (a.y == b.y) && (a.z == b.z) && (a.w == b.w);
}


template< typename VectorType >
inline typename VectorType::real_t distance(const VectorType &a, const VectorType &b)
{
    return magnitude( a - b );
}

template< typename VectorType >
inline typename VectorType::real_t distance_sq(const VectorType &a, const VectorType &b)
{
    const VectorType c(a - b);
    return dot( c, c );
}

template< typename VectorType>
inline VectorType normalize(const VectorType &a)
{
    return a / magnitude(a);
}

