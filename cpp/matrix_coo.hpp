// matrix_coo.hpp
//
// James Folberth
// Spring 2015

#ifndef _MATRIX_COO_HPP_
#define _MATRIX_COO_HPP_

#include <algorithm>
#include <cassert>
#include <iostream>
#include <iomanip>
#include <vector>

#include "utils.hpp"
#include "matrix_crs.hpp"

using namespace std;

// Forward declare to make function declaration possible
template<typename T>
class matrix_coo;

template<typename T>
class matrix_crs;

template<typename T>
ostream& operator<<(ostream& os, const matrix_coo<T>& mat); 

template <typename T>
class matrix_coo {
   public:
      using value_type = T;

      vector<unsigned> row_ind;
      vector<unsigned> col_ind;
      vector<T> val;
      size_t m;          // number of rows
      size_t n;          // number of cols

      //////////////////////////////////
      // Construction and Destruction //
      //////////////////////////////////
      matrix_coo<T>() = default;

      matrix_coo<T>(vector<unsigned>& init_row_ind,
                 vector<unsigned>& init_col_ind,
                 vector<T>& init_val,
                 size_t init_m=0, size_t init_n=0);

      matrix_coo<T>(const matrix_coo<T>& ) = default;
      matrix_coo<T>& operator=(const matrix_coo<T>& ) = default;

      ~matrix_coo<T>() = default;

      /////////////////////
      // Type conversion //
      /////////////////////
      matrix_crs<T> to_crs(void);


      ////////////
      // Output //
      ////////////
      // Defined below, not in matrix_coo.cpp
      friend ostream& operator<< <>(ostream& os, const matrix_coo& mat); 

      void print_full(void);
 
   private:
      void sort_inds(void);
 
};

template <typename T>
ostream& operator<<(ostream& os, const matrix_coo<T>& mat) {

   os << "Coordinate Sparse (COO) (rows = " << mat.m << ", cols = " << 
      mat.n << ", nnz = " << mat.val.size() << ")" << endl; // << endl;

   if ( _DEBUG_ >=1 ) {
      os << "debug 1: COO ostream printing" << endl;
      os << "row_ind: ";
      for (unsigned i=0; i < mat.row_ind.size(); ++i)
         cout << mat.row_ind[i] << "  ";
      os << endl;
      os << "col_ind: ";
      for (unsigned i=0; i < mat.col_ind.size(); ++i)
         cout << mat.col_ind[i] << "  ";
      os << endl;
      os << "val:     ";
      for (unsigned i=0; i < mat.val.size(); ++i)
         cout << setprecision(_PRINT_SPARSE_PREC_) << mat.val[i] << "  ";
      os << endl;
   }

   for ( unsigned i=0; i < mat.row_ind.size(); ++i) {
      os << "  (" << mat.row_ind[i] << ", " << mat.col_ind[i]
         << ") -> " << setprecision(_PRINT_SPARSE_PREC_)
         << mat.val[i] << endl;
   }

   return os;
}

#endif
