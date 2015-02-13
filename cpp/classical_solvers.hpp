// classical_solvers.hpp
//
// James Folberth
// Spring 2015

#ifndef _CLASSICAL_SOLVERS_HPP_
#define _CLASSICAL_SOLVERS_HPP_

#include <valarray>

#include "matrix_crs.hpp"
#include "utils.hpp"

#define _WJ_DEFAULT_W_ 2./3.

#define _WJ_DEFAULT_RESID_TOL_ 10e-8
#define _GS_DEFAULT_RESID_TOL_ 10e-8
#define _RBGS_DEFAULT_RESID_TOL_ 10e-8

#define _WJ_MAX_ITR_ 100
#define _GS_MAX_ITR_ 100
#define _RBGS_MAX_ITR_ 100

using namespace std;

////////
// WJ //
////////
template<typename T>
valarray<T> wjacobi(const matrix_crs<T>& A, const valarray<T>& f, 
      const valarray<T>& v0, const T w=_WJ_DEFAULT_W_,
      const T resid_tol=_WJ_DEFAULT_RESID_TOL_, int& num_itr=-1);

// specify resid (or default)
template<typename T>
valarray<T> wjacobi(const matrix_crs<T>& A, const valarray<T>& f, 
      const T w=_WJ_DEFAULT_W_, const T resid_tol=_WJ_DEFAULT_RESID_TOL_,
      int& num_itr=-1);

template<typename T>
valarray<T> wjacobi(const matrix_crs<T>& A, const valarray<T>& f, 
      const T w=_WJ_DEFAULT_W_, const T resid_tol=_WJ_DEFAULT_RESID_TOL_);

// specify num_itr (or default)
template<typename T>
valarray<T> wjacobi(const matrix_crs<T>& A, const valarray<T>& f, 
      const valarray<T>& v0, const T w=_WJ_DEFAULT_W_, int& num_itr = -1);

template<typename T>
valarray<T> wjacobi(const matrix_crs<T>& A, const valarray<T>& f, 
      const T w=_WJ_DEFAULT_W_, int& num_itr = -1);

template<typename T>
void wjacobi_it(const matrix_crs<T>& A, const valarray<T>& f,
      const valarray<T>& v0, valarray<T>& v1, const T w);


////////
// GS //
////////
// main driver
template<typename T>
valarray<T> gauss_seidel(const matrix_crs<T>& A, const valarray<T>& f, 
      const valarray<T>& v0, const T resid_tol=_GS_DEFAULT_RESID_TOL_,
      int& num_itr = -1);

// specify resid (or default)
template<typename T>
valarray<T> gauss_seidel(const matrix_crs<T>& A, const valarray<T>& f,
      const T resid_tol=_GS_DEFAULT_RESID_TOL_, int& num_itr=-1);

template<typename T>
valarray<T> gauss_seidel(const matrix_crs<T>& A, const valarray<T>& f,
      const T resid_tol=_GS_DEFAULT_RESID_TOL_);

// specify num_itr (or default)
template<typename T>
valarray<T> gauss_seidel(const matrix_crs<T>& A, const valarray<T>& f, 
      const valarray<T>& v0, int& num_itr = -1);

template<typename T>
valarray<T> gauss_seidel(const matrix_crs<T>& A, const valarray<T>& f, 
      int& num_itr = -1);

template<typename T>
void gauss_seidel_it(const matrix_crs<T>& A, const valarray<T>& f,
      valarray<T>& v);


//////////
// RBGS //
//////////
// main driver
template<typename T>
valarray<T> rbgauss_seidel(const matrix_crs<T>& A, const valarray<T>& f, 
      const valarray<T>& v0, const T resid_tol=_GS_DEFAULT_RESID_TOL_,
      int& num_itr = -1);

// specify resid (or default)
template<typename T>
valarray<T> rbgauss_seidel(const matrix_crs<T>& A, const valarray<T>& f,
      const T resid_tol=_GS_DEFAULT_RESID_TOL_, int& num_itr=-1);

template<typename T>
valarray<T> rbgauss_seidel(const matrix_crs<T>& A, const valarray<T>& f,
      const T resid_tol=_GS_DEFAULT_RESID_TOL_);

// specify num_itr (or default)
template<typename T>
valarray<T> rbgauss_seidel(const matrix_crs<T>& A, const valarray<T>& f, 
      const valarray<T>& v0, int& num_itr = -1);

template<typename T>
valarray<T> rbgauss_seidel(const matrix_crs<T>& A, const valarray<T>& f, 
      int& num_itr = -1);

template<typename T>
void rbgauss_seidel_it(const matrix_crs<T>& A, const valarray<T>& f,
      valarray<T>& v);


#endif
