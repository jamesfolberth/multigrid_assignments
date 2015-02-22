// multigrid.hpp
//
// James Folberth
// Spring 2015

#ifndef _MULTIGRID_HPP_
#define _MULTIGRID_HPP_

#include <functional> // for passing capture-full lambda functions
#include <vector>
#include <valarray>

#include "matrix_crs.cpp"
#include "utils.hpp"
#include "classical_solvers.hpp"

using namespace std;

///////////
// Level //
///////////
// {{{

template<typename T>
class level;

template<typename T>
class level {
   public:
      using value_type = T;

      matrix_crs<T> A;
      valarray<T> v;
      valarray<T> f;
      
      matrix_crs<T> P; // prolongation operator
      matrix_crs<T> R; // restriction operator

      // smoother_ip
      function<void(const matrix_crs<T>&, const valarray<T>&,
            valarray<T>&, unsigned)> smoother_ip;

      // construction
      level<T>() = default;

      // destruction, move, copy
      level<T>(const level<T>& ) = default;
      level<T>& operator=(const level<T>& ) = default;
      level<T>(level<T>&& ) = default;
      level<T>& operator=(level<T>&& ) = default;
      ~level<T>() = default;

};

template<typename T>
vector<level<T>> build_levels(function<matrix_crs<T>(unsigned)> build_A,
      valarray<T> f, function<matrix_crs<T>(unsigned)> build_P, 
      function<matrix_crs<T>(unsigned)> build_R,
      function<void(const matrix_crs<T>&, const valarray<T>&,
         valarray<T>&, unsigned)> smoother_ip,
      unsigned L, unsigned Lx, const valarray<T>& v0);

template<typename T>
vector<level<T>> build_levels(function<matrix_crs<T>(unsigned)> build_A,
      valarray<T> f, function<matrix_crs<T>(unsigned)> build_P, 
      function<matrix_crs<T>(unsigned)> build_R,
      function<void(const matrix_crs<T>&, const valarray<T>&,
         valarray<T>&, unsigned)> smoother_ip,
      unsigned L, unsigned Lx);

// }}}

////////////
// Cycles //
////////////
// {{{
template<typename T>
void mucycle(vector<level<T>>& levels, typename vector<level<T>>::iterator it,
      unsigned nu1, unsigned nu2, unsigned mu=1);

template<typename T>
void vcycle(vector<level<T>>& levels, typename vector<level<T>>::iterator it,
      unsigned nu1, unsigned nu2);

// }}}

/////////////////////////
// Intergrid Operators //
/////////////////////////
// {{{
template<typename T>
matrix_crs<T> operator_1d_interp_lin(unsigned lm1);

template<typename T>
matrix_crs<T> operator_1d_restrict_inj(unsigned l);

template<typename T>
matrix_crs<T> operator_1d_restrict_full(unsigned l);

// }}}


#endif
