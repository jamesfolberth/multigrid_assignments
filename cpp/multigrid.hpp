// multigrid.hpp
//
// James Folberth
// Spring 2015

#ifndef _MULTIGRID_HPP_
#define _MULTIGRID_HPP_

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
