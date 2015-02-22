// multigrid.cpp
//
// James Folberth
// Spring 2015

#include "multigrid.hpp"

using namespace std;

///////////
// Level //
///////////


/////////////////////////
// Intergrid Operators //
/////////////////////////
// {{{
template<typename T>
matrix_crs<T> operator_1d_interp_lin(unsigned lm1) {
// form the linear operator P that performs linear interpolation from grid
// lm1 = l-1 to grid l.  Number of points on the initial grid is 2^(lm1)-1

   unsigned m = pow(2,lm1+1)-1, n = pow(2,lm1)-1;

   // Build vectors for the i,j,value tuples
   // There are three entries per column; there are n columns
   vector<unsigned> row_ind(3*n), col_ind(3*n);
   vector<T> val(3*n);

   for (unsigned col=0; col < n; ++col) {
      row_ind[3*col]   = 2*col;
      row_ind[3*col+1] = 2*col+1;
      row_ind[3*col+2] = 2*col+2;

      col_ind[3*col]   = col;
      col_ind[3*col+1] = col;
      col_ind[3*col+2] = col;

      val[3*col]   = static_cast<T>(0.5);
      val[3*col+1] = static_cast<T>(1);
      val[3*col+2] = static_cast<T>(0.5);
   }
 
   // build CRS matrix and return
   return matrix_crs<T>(row_ind, col_ind, val, m, n);
}

template<typename T>
matrix_crs<T> operator_1d_restrict_inj(unsigned l) {
// form the linear operator R that restricts from grid l to grid l-1 using
// injection.  Number of grid points on the initial grid is 2^l-1

   unsigned m = pow(2,l-1)-1, n = pow(2,l)-1;

   // Build vectors for the i,j,value tuples
   // There is only one element per row of the matrix
   vector<unsigned> row_ind(m), col_ind(m);
   vector<T> val(m);

   for (unsigned row=0; row < m; ++row) {
      row_ind[row] = row;
      col_ind[row] = 2*row+1;
      val[row] = static_cast<T>(1.);
   }

   // build CRS matrix and return
   return matrix_crs<T>(row_ind, col_ind, val, m, n);
}

template<typename T>
matrix_crs<T> operator_1d_restrict_full(unsigned l) {
// form the linear operator R that restricts from grid l to grid l-1 using
// full weighting.  Number of grid points on the initial grid is 2^l-1

   unsigned m = pow(2,l-1)-1, n = pow(2,l)-1;

   // Build vectors for the i,j,value tuples
   // There are three elements per row of the matrix
   vector<unsigned> row_ind(3*m), col_ind(3*m);
   vector<T> val(3*m);

   for (unsigned row=0; row < m; ++row) {
      row_ind[3*row] = row;
      row_ind[3*row+1] = row;
      row_ind[3*row+2] = row;

      col_ind[3*row] = 2*row;
      col_ind[3*row+1] = 2*row+1;
      col_ind[3*row+2] = 2*row+2;

      val[3*row] = static_cast<T>(0.25);
      val[3*row+1] = static_cast<T>(0.5);
      val[3*row+2] = static_cast<T>(0.25);
   }

   // build CRS matrix and return
   return matrix_crs<T>(row_ind, col_ind, val, m, n);
}
// }}}


// Force instatiation
// double
template matrix_crs<double> operator_1d_interp_lin<double>(unsigned);
template matrix_crs<double> operator_1d_restrict_inj<double>(unsigned);
template matrix_crs<double> operator_1d_restrict_full<double>(unsigned);


