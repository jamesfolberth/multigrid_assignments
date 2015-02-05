// matrix_coo.cpp
//
// James Folberth
// Spring 2015

#include "matrix_coo.hpp"

using namespace std;

//////////////////////////////////
// Construction and destruction //
//////////////////////////////////
template<typename T>
matrix_coo<T>::matrix_coo(
      vector<unsigned>& init_row_ind,
      vector<unsigned>& init_col_ind,
      vector<T>& init_val,
      size_t init_m, size_t init_n) {

   row_ind = init_row_ind;
   col_ind = init_col_ind;
   val = init_val;

   if (val.size() == 0) {
      if (init_m == 0 && init_n == 0) {
         cerr << "error: matrix_coo: Can't construct an empty (zero) matrix "
              << "without specifying m and n!" << endl;
         exit(-1);
      }
      m = init_m; n = init_n;
      return;
   }

   unsigned max_row = 1+*max_element(row_ind.begin(), row_ind.end());
   unsigned max_col = 1+*max_element(col_ind.begin(), col_ind.end());
   m = (init_m < max_row) ? max_row : init_m;
   n = (init_n < max_col) ? max_col : init_n;

   assert(row_ind.size() == col_ind.size());
   assert(col_ind.size() == val.size());

   this->clean();
}

// TODO should accept a combine function (default to lambda add) like
// Julia's CSC
template<typename T>
void matrix_coo<T>::clean(void) {
// Reorder the row/col/vals so that they are in row-major order

   // Copy row inds, col inds, vals to vector of (i,j,val) tuples
   // Then sort vector of tuples with stdlib sort
   // Then overwrite values and remove duplicates
   // I don't know how fast this is
   
   struct sort_tuple {
      unsigned i;
      unsigned j;
      T val;
   };

   vector<sort_tuple> sort_me;
   sort_me.resize(row_ind.size());
   for (unsigned i=0; i<sort_me.size(); ++i) {
      sort_me[i].i = row_ind[i];
      sort_me[i].j = col_ind[i];
      sort_me[i].val = val[i];
   }

   sort(sort_me.begin(), sort_me.end(), 
         [&] (const sort_tuple& lhs, const sort_tuple& rhs) -> bool {
            if (lhs.i < rhs.i) {
               return 1;
            }
            else if (lhs.i == rhs.i) {
               if (lhs.j < rhs.j) return 1;
               else return 0; // ignore repeats
            }
            else return 0;
         });

   // assign to class member and deal with duplicate entries by adding
   unsigned old_row=-1, old_col=-1, num_dups=0;
   for (unsigned i=0; i<sort_me.size(); ++i) {

      if ( abs(sort_me[i].val) < _ELEMENT_ZERO_TOL_ ) {
         num_dups += 1;
         continue;
      }

      // if not a duplicate (case i=0 always passes, since it's first)
      if (old_row != sort_me[i].i || old_col != sort_me[i].j) {

         row_ind[i-num_dups] = sort_me[i].i;
         col_ind[i-num_dups] = sort_me[i].j;
         val[i-num_dups] = sort_me[i].val;

         old_row = sort_me[i].i;
         old_col = sort_me[i].j;
      }
      
      else {
         num_dups += 1;
         val[i-num_dups] += sort_me[i].val; // TODO combine function
         continue;
      }
   }

   if (num_dups > 0) {
      cerr << "warning: matrix_coo:clean: duplicate entries found or "
           << "zeros; combining entries by adding and removing zero entries"
           << endl;
   }

   row_ind.resize(row_ind.size()-num_dups);
   col_ind.resize(col_ind.size()-num_dups);
   val.resize(val.size()-num_dups);
}



// TODO maybe move to seperate file?
// special forms
template<typename T>
matrix_coo<T> eye_coo(unsigned m, unsigned n) {
   vector<unsigned> rind(MIN(m,n)), cind(MIN(m,n));
   vector<T> val(MIN(m,n), static_cast<T>(1)); // fill val with 1s on construct

   for (unsigned i=0; i < rind.size(); ++i) {
      rind[i] = i;
      cind[i] = i;
   }

   return matrix_coo<T>(rind, cind, val, m, n);
}


////////////////
// Operations //
////////////////
// scalar
template<typename T>
matrix_coo<T>& matrix_coo<T>::operator+=(const T& value) {
   for (auto it=val.begin(); it != val.end(); ++it) {
      *it += value;
   }
   return *this;
}

template<typename T>
matrix_coo<T>& matrix_coo<T>::operator-=(const T& value) {
   for (auto it=val.begin(); it != val.end(); ++it) {
      *it -= value;
   }
   return *this;
}

template<typename T>
matrix_coo<T>& matrix_coo<T>::operator*=(const T& value) {
   for (auto it=val.begin(); it != val.end(); ++it) {
      *it *= value;
   }
   return *this;
}

template<typename T>
matrix_coo<T>& matrix_coo<T>::operator/=(const T& value) {
   for (auto it=val.begin(); it != val.end(); ++it) {
      *it /= value;
   }
   return *this;
}

/////////////////////
// Type conversion //
/////////////////////
template<typename T>
matrix_crs<T> matrix_coo<T>::to_crs(void) {
   return matrix_crs<T>(row_ind, col_ind, val, m, n);
}


////////////
// Output //
////////////
template<typename T>
void matrix_coo<T>::print_full(void) {
   this->to_crs().print_full();
}


// Some template functions are defined in the header
// I couldn't get things to work any other way

// Force instantiation for specific types
// double
template class matrix_coo<double>;
template matrix_coo<double> eye_coo<double>(unsigned, unsigned);

