#ifndef __VECTOR_H
#define __VECTOR_H

#include <iostream>

template <class F>
struct Vector {
  F x, y, z;
};


template <class F>
inline __device__ __host__ Vector<F>
operator+(const Vector<F> &v1,
	  const Vector<F> &v2)
{
  Vector<F> r;

  r.x = v1.x + v2.x;
  r.y = v1.y + v2.y;
  r.z = v1.z + v2.z;

  return r;
}

template <class F>
inline __device__ __host__ Vector<F>
operator-(const Vector<F> &v1,
	  const Vector<F> &v2)
{
  Vector<F> r;

  r.x = v1.x - v2.x;
  r.y = v1.y - v2.y;
  r.z = v1.z - v2.z;

  return r;
}

template <class F>
inline __device__ __host__ F
operator*(const Vector<F> &v1,
	  const Vector<F> &v2)
{
  return v1.x*v2.x + v1.y*v2.y + v1.z*v2.z;
}

template <class F>
inline __device__ __host__ Vector<F>
operator*(F a, const Vector<F> &v)
{
  Vector<F> r;

  r.x = a * v.x;
  r.y = a * v.y;
  r.z = a * v.z;

  return r;
}

template <class F>
inline __device__ __host__ Vector<F>
operator&(const Vector<F> &v1,
	  const Vector<F> &v2)
{
  Vector<F> r;

  r.x = v1.y*v2.z - v1.z*v2.y;
  r.y = v1.z*v2.x - v1.x*v2.z;
  r.z = v1.x*v2.y - v1.y*v2.x;

  return r;
}

template <class F>
inline __device__ __host__ F
norm(const Vector<F> &v)
{
  return sqrt(v * v);
}

template <class F>
inline __device__ __host__ Vector<F>
normalize(const Vector<F> &v)
{
  return 1/norm(v) * v;
}

template <class F>
inline std::ostream & operator<<(std::ostream &o, const Vector<F> &v)
{
  o << "x: " << v.x << " y: " << v.y << " z: " << v.z ;
  return o;
}

#endif

