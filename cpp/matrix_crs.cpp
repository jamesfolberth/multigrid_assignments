// matrix_crs.cpp
//
// James Folberth
// Spring 2015

#include "matrix_crs.hpp"

using namespace std;

//////////////////////////////////
// Construction and Destruction //
//////////////////////////////////
template<typename T>
matrix_crs<T>::matrix_crs(
      vector<unsigned>& init_row_ind,
      vector<unsigned>& init_col_ind,
      vector<T>& init_val,
      size_t init_m, size_t init_n) {
// Construct a CRS matrix from vectors of row index, column index, and values
// We first build a COO matrix to sort and deal with duplicate entries, and
// also determine the size of the matrix.
// Then we build the vector of row pointers, taking into account zero rows

   // Build COO matrix (which will sort indexes)
   matrix_coo<T> coomat(init_row_ind, init_col_ind, init_val, init_m, init_n); 

   if ( _DEBUG_ >= 2) {
      cout << "debug 2: CRS matrix construction" << endl;
      cout << coomat << endl;
   }

   col_ind = coomat.col_ind;
   val = coomat.val;
   m = coomat.m;
   n = coomat.n;

   // Build vector of row pointers
   row_ptr.resize(m+1);
   row_ptr[m] = val.size();
   unsigned row_ptr_ind = 0;

   // Test for matrix of zeros
   if (val.size() == 0) {
      fill(row_ptr.begin(), row_ptr.end(), 0);
      return;
   }

   // Handle initial rows of zeros
   for (unsigned j=0; j<coomat.row_ind[0]; ++j)
      row_ptr[row_ptr_ind + j] = 0;
   row_ptr_ind = coomat.row_ind[0];

   for (unsigned i=1; i<col_ind.size(); ++i) {

      // nothing to do; keep indexing along row
      if (coomat.row_ind[i-1] == coomat.row_ind[i]) {

         // unless it's the last row
         if (coomat.row_ind[i] == m) {
            row_ptr[row_ptr_ind+1] = i;
         }
         continue;
      }

      // jumping down to next column
      else if (coomat.row_ind[i-1] == coomat.row_ind[i] - 1) {
         row_ptr[row_ptr_ind+1] = i;
         row_ptr_ind += 1;
      }

      // we have a row or rows of zeros
      // Repeat last index of col_ind for each row we skip
      // Also put in the next row_ptr for the next non-zero row
      else if (coomat.row_ind[i-1] < coomat.row_ind[i] - 1 ) {
         for (unsigned j=1; j<coomat.row_ind[i]-row_ptr_ind; ++j) {
            row_ptr[row_ptr_ind + j] = i;
         }
         row_ptr_ind = coomat.row_ind[i];
         row_ptr[row_ptr_ind] = i;
      }

      else {
         cerr << "error: matrix_crs.cpp::matrix_crs: error constructing "
              << "CRS matrix" << endl;
         exit(-1);
      }
   }

   // Handle rows of zeros at the end of the matrix
   for (unsigned r=*(coomat.row_ind.end()-1)+1; r<m; ++r) {
      row_ptr[r] = col_ind.size();
   }
}


/////////////////////
// Type conversion //
/////////////////////
template<typename T>
matrix_coo<T> matrix_crs<T>::to_coo(void) {
   vector<unsigned> new_row_ind(col_ind.size());

   for (unsigned i=0; i<row_ptr.size()-1; ++i) {
      fill(new_row_ind.begin()+row_ptr[i],
           new_row_ind.begin()+row_ptr[i+1],
           i);
   }

   return matrix_coo<T>(new_row_ind, col_ind, val, m, n);
}


template<typename T>
matrix_crs<T> matrix_crs<T>::to_crs(void) {
   return *this;
}



////////////
// Output //
////////////
template<typename T>
void matrix_crs<T>::print_full(void) {
// Print sparse matrix as a dense array

   unsigned ind=0;
   for (unsigned r=0; r<m; ++r) {
      for (unsigned c=0; c<n; ++c) {

         // we've hit a nonzero in this row
         if (col_ind.size() > 0 && c == col_ind[ind] && ind < row_ptr[r+1] ) { 
            cout << ' ' << setw(_PRINT_FULL_WIDTH_) << setfill(' ')
                 << setprecision(_PRINT_FULL_PREC_)
                 << static_cast<double>(val[ind]);
            ind += 1;
         }

         // print a zero
         else {
            cout << ' ' << setw(_PRINT_FULL_WIDTH_) << setfill(' ')
                 << setprecision(_PRINT_FULL_PREC_)
                 << 0e0;
         }

         // end of row, so print endl
         if (c == n-1) cout << endl;
      }
   }

}


// Some template functions are defined in the header
// I couldn't get things to work any other way

// Force instantiation for specific types
//template class matrix_crs<float>;
template class matrix_crs<double>;

