// matrix_crs.hpp
//
// James Folberth
// Spring 2015

#ifndef _MATRIX_CRS_HPP_
#define _MATRIX_CRS_HPP_

#include <algorithm>
#include <cassert>
#include <iostream>
#include <vector>

#include "utils.hpp"
#include "matrix_coo.hpp"

using namespace std;

// Forward declare to make function declaration possible
template<typename T>
class matrix_crs;

template<typename T>
ostream& operator<<(ostream& os, const matrix_crs<T>& mat); 

template <typename T>
class matrix_crs {
   public:
      using value_type = T;

      vector<unsigned> row_ptr;
      vector<unsigned> col_ind;
      vector<T> val;
      size_t m;          // number of rows
      size_t n;          // number of cols

      // Construction and Destruction
      matrix_crs<T>() = default;

      matrix_crs<T>(vector<unsigned>& init_row_ind,
                 vector<unsigned>& init_col_ind,
                 vector<T>& init_val,
                 size_t init_m=0, size_t init_n=0);

      ~matrix_crs<T>() = default;// {cout << "delete me!" << endl;};

      // Type conversion
      //matrix_coo<T>& to_coo();

      // Output
      // Defined below, not in matrix_crs.cpp
      friend ostream& operator<< <>(ostream& os, const matrix_crs& mat); 

      void print_full(void);
 
   private:
      void sort_inds(void);
 
};

template <typename T>
ostream& operator<<(ostream& os, const matrix_crs<T>& mat) {

   os << "Compressed Row Storage (CRS) (rows = " << mat.m << ", cols = " << 
      mat.n << ", nnz = " << mat.val.size() << ")" << endl; // << endl;

   if ( _DEBUG_ >= 1 ) {
      os << "row_ptr: ";
      for (unsigned i=0; i<mat.row_ptr.size(); ++i)
         os << mat.row_ptr[i] << "  ";
      os << endl;
      os << "col_ind: ";
      for (unsigned i=0; i<mat.col_ind.size(); ++i)
         os << mat.col_ind[i] << "  ";
      os << endl;
      os << "val:     ";
      for (unsigned i=0; i<mat.val.size(); ++i)
         os << mat.val[i] << "  ";
      os << endl;
   }

   for (unsigned i=0; i<mat.m; ++i) {
      for (unsigned j=mat.row_ptr[i]; j<mat.row_ptr[i+1]; ++j) {
         cout << "  (" << i << ", " << mat.col_ind[j] << ") -> " 
            << setprecision(_PRINT_SPARSE_PREC_) << mat.val[j] << endl;
      }
   }

   // TODO implement me
   return os;
}

#endif
