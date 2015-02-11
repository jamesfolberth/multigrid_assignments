// classical_solvers.hpp
//
// James Folberth
// Spring 2015

#ifndef _CLASSICAL_SOLVERS_HPP_
#define _CLASSICAL_SOLVERS_HPP_

#include <valarray>

#include "matrix_crs.hpp"
#include "utils.hpp"

using namespace std;

template<typename T>
valarray<T> wjacobi(const matrix_crs<T>& A, const valarray<T>& f, 
                    const T w=2./3., int num_itr = -1);

template<typename T>
void wjacobi_it(const matrix_crs<T>& A, const valarray<T>& f,
                       const valarray<T>& v0, valarray<T>& v1,
                       const T w);


#endif
