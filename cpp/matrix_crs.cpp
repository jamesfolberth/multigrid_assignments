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


// TODO do this in place on the CRS data instead of converting and cleaning
template<typename T>
void matrix_crs<T>::clean(void) {
   matrix_coo<T> temp = this->to_coo();
   *this = temp.to_crs();
}


// TODO maybe move to seperate file?
// special forms
template<typename T>
matrix_crs<T> eye_crs(unsigned m, unsigned n) {
   vector<unsigned> rind(MIN(m,n)), cind(MIN(m,n));
   vector<T> val(MIN(m,n), static_cast<T>(1)); // fill val with 1s on construct

   for (unsigned i=0; i < rind.size(); ++i) {
      rind[i] = i;
      cind[i] = i;
   }

   return matrix_crs<T>(rind, cind, val, m, n);
}


////////////////
// Operations //
////////////////
// scalar
template<typename T>
matrix_crs<T>& matrix_crs<T>::operator+=(const T& value) {
   for (auto it=val.begin(); it != val.end(); ++it) {
      *it += value;
   }
   return *this;
}

template<typename T>
matrix_crs<T>& matrix_crs<T>::operator-=(const T& value) {
   for (auto it=val.begin(); it != val.end(); ++it) {
      *it -= value;
   }
   return *this;
}

template<typename T>
matrix_crs<T>& matrix_crs<T>::operator*=(const T& value) {
   for (auto it=val.begin(); it != val.end(); ++it) {
      *it *= value;
   }
   return *this;
}

template<typename T>
matrix_crs<T>& matrix_crs<T>::operator/=(const T& value) {
   for (auto it=val.begin(); it != val.end(); ++it) {
      *it /= value;
   }
   return *this;
}

// matrix add
// *this += B;
template<typename T>
matrix_crs<T>& matrix_crs<T>::operator+=(const matrix_crs<T>& B) {
   unsigned this_col_ind;

   // check sizes
   if ( m != B.m || n != B.n ) {
      cerr << "error: matrix_crs:+=: matrix sizes do not match - ("
           << m << "," << n << ") vs. ("
           << B.m << "," << B.n << ")" << endl;
      exit(-1);
   }

   // if empty mat, do nothing
   if ( val.size() == 0 ) {
      return *this;
   }

   // loop through rows
   for (unsigned row=0; row < m; ++row) {
      this_col_ind = row_ptr[row];
    
      // for each row, loop through column indexes of B
      // keep a column index of *this to compare with the column index of B
      for (unsigned B_col_ind = B.row_ptr[row];
            B_col_ind < B.row_ptr[row+1];
            ++B_col_ind) {

         // increase this_col_ind until we hit the end of the row
         // or an entry of B 
         while ( this_col_ind < row_ptr[row+1] && 
                 col_ind[this_col_ind] < B.col_ind[B_col_ind] )
            ++this_col_ind; 

         // this_col_ind still in the row
         if ( this_col_ind < row_ptr[row+1] ) {

            // entry in *this and B - add them
            if ( col_ind[this_col_ind] == B.col_ind[B_col_ind] ) {
               val[this_col_ind] += B.val[B_col_ind];

               // if val below _ELEMENT_ZERO_TOL_, remove it
               if ( abs(val[this_col_ind]) < _ELEMENT_ZERO_TOL_ ) {
                  col_ind.erase(col_ind.begin()+this_col_ind);
                  val.erase(val.begin()+this_col_ind);

                  for (unsigned r=row+1; r<=m; ++r) {
                     row_ptr[r] -= 1;
                  }
               }
               continue;
            }
            
            // entry in B before *this - insert entry 
            // Note: this could be an else, instead of else if; the cases
            // are exhaustive
            else if ( col_ind[this_col_ind] > B.col_ind[B_col_ind] ) {
               col_ind.insert(col_ind.begin()+this_col_ind,
                     B.col_ind[B_col_ind]);
               val.insert(val.begin()+this_col_ind,
                     B.val[B_col_ind]);

               // fix row pointers
               for (unsigned r=row; r <= m; ++r) {
                  row_ptr[r] += 1;
               }
               continue;
            }

            else {
               cerr << "error: matrix_crs:+=: I have reached an impossible "
                    << "state!  Dying!" << endl;
               exit(-1);
            }
         }

         // this_col_ind is not still in the row (i.e. it's in the next row)
         // still need to insert entry of B. 
         else {
            // can't have entry in *this and B, since past row of *this
            // Thus, insert entry of B
            col_ind.insert(col_ind.begin()+this_col_ind,
                  B.col_ind[B_col_ind]);
            val.insert(val.begin()+this_col_ind,
                  B.val[B_col_ind]);

            // fix row pointers
            for (unsigned r=row; r <= m; ++r) {
               row_ptr[r] += 1;
            }
         }
      }
   }

   return *this;
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


//template<typename T>
//matrix_crs<T> matrix_crs<T>::to_crs(void) {
//   return *this;
//}



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

   cout << endl;
}


// Some template functions are defined in the header
// I couldn't get things to work any other way

// Force instantiation for specific types
// double
template class matrix_crs<double>;
template matrix_crs<double> eye_crs<double>(unsigned,unsigned);

