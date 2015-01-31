// matrix_coo.cpp
//
// James Folberth
// Spring 2015

#include "matrix_coo.hpp"

using namespace std;

// Construction and destruction
template<typename T>
matrix_coo<T>::matrix_coo(
      vector<unsigned>& init_row_ind,
      vector<unsigned>& init_col_ind,
      vector<T>& init_val,
      size_t init_m, size_t init_n) {

   row_ind = init_row_ind;
   col_ind = init_col_ind;
   val = init_val;

   unsigned max_row = 1+*max_element(row_ind.begin(), row_ind.end());
   unsigned max_col = 1+*max_element(col_ind.begin(), col_ind.end());
   m = (init_m < max_row) ? max_row : init_m;
   n = (init_n < max_col) ? max_col : init_n;

   assert(row_ind.size() == col_ind.size());
   assert(col_ind.size() == val.size());

   this->sort_inds();
}


template<typename T>
void matrix_coo<T>::sort_inds(void) {
// Reorder the row/col/vals so that they are in row-major order

   // Copy row inds, col inds, vals to vector of (i,j,val) tuples
   // Then sort vector of tuples with stdlib sort
   // Then overwrite values
   // This is going to be slow

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

   for (unsigned i=0; i<sort_me.size(); ++i) {
      row_ind[i] = sort_me[i].i;
      col_ind[i] = sort_me[i].j;
      val[i] = sort_me[i].val;
   }
}


// Some template functions are defined in the header
// I couldn't get things to work any other way

// Force instantiation for specific types
//template class matrix_coo<float>;
template class matrix_coo<double>;

