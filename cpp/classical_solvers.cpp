// classical_solvers.cpp
//
// James Folberth
// Spring 2015

#include "classical_solvers.hpp"

using namespace std;

/////////////////////
// Weighted Jacobi //
/////////////////////
// {{{
template<typename T>
valarray<T> wjacobi(const matrix_crs<T>& A, const valarray<T>& f, const T w,
      int num_itr) {
   // weighted Jacobi method driver
   // pass in the matrix, RHS vector
   // weight is optional.  defaults to 2./3.
   // num_itr is the number of iterations to do.  defaults to -1, which will
   //   run until we reach the maximum number of iterations or reach tolerance
   //   On return, num_itr will have the number of iterations completed.

   unsigned max_itr = 100;
   T resid_tol = static_cast<T>(10e-8);
   
   unsigned num_its_done, num_its_todo;
   
   // check sizes
   if ( A.m != f.size() ) {
      cerr << "error: classical_solvers:wjacobi: dimension mismatch" << endl;
      exit(-1);
   }

   // random initial guess
   valarray<T> v = rand_vec<T>(A.n,-1.,1.);
   valarray<T> v_prev = v;
   valarray<T> resid = f-A*v;

   // number of iterations to do
   if ( num_itr == -1 ) { // just do normal solve until we reach tol or max_itr
   
      // do weighted Jacobi iterations
      // norm(...,0) is infinity norm
      num_its_done = 0;
      while ( norm(resid,0) > resid_tol && num_its_done < max_itr ) {
         v_prev = v;
         
         wjacobi_it(A,f,v_prev,v,w);
   
         resid = f-A*v;
         ++num_its_done;
      }
   }

   else if ( num_itr >= 0 ) { // do exactly num_itr iterations
      // do weighted Jacobi iterations
      // norm(...,0) is infinity norm
      num_its_todo = static_cast<unsigned>(num_itr);
      num_its_done = 0;
      while ( num_its_done < num_its_todo ) {
         v_prev = v;
         
         wjacobi_it(A,f,v_prev,v,w);
   
         ++num_its_done;
      }
   }
   
   else {
      cerr << "error: classical_solvers:wjacobi: bad input number of "
           << "iterations" << endl;
      exit(-1);
   }

   num_itr = num_its_done;
   return v;
}

template<typename T>
void wjacobi_it(const matrix_crs<T>& A, const valarray<T>& f,
             const valarray<T>& v0, valarray<T>& v1, const T w) {
   // Do one iteration of weighted Jacobi
   T LpUv, ajj=0.;

   // sweep down rows of A
   for (size_t row=0; row < A.m; ++row) {
      LpUv = 0.;
   
      // sweep across column.  accumulated (L+U)*v0 and find a_{jj}
      for (size_t ptr=A.row_ptr[row]; ptr < A.row_ptr[row+1]; ++ptr) {
         if ( A.col_ind[ptr] != row ) {// (L+U)*v
            LpUv -= A.val[ptr] * v0[A.col_ind[ptr]];
         }
         else {// get a_{jj}
            ajj = A.val[ptr];
         }
      }
   
      // assemble for the update
      v1[row] = (1.-w)*v0[row]+w*(LpUv+f[row])/ajj;
   }
}

// }}}

//////////////////
// Gauss-Seidel //
//////////////////
// {{{
template<typename T>
valarray<T> gauss_seidel(const matrix_crs<T>& A, const valarray<T>& f,
                         int num_itr) {
   // Gauss Seidel method driver
   // pass in the matrix, RHS vector
   // weight is optional.  defaults to 2./3.
   // num_itr is the number of iterations to do.  defaults to -1, which will
   //   run until we reach the maximum number of iterations or reach tolerance
   //   On return, num_itr will have the number of iterations completed.

   unsigned max_itr = 100;
   T resid_tol = static_cast<T>(10e-8);
   
   unsigned num_its_done, num_its_todo;
   
   // check sizes
   if ( A.m != f.size() ) {
      cerr << "error: classical_solvers:gauss_seidel: dimension mismatch"
           << endl;
      exit(-1);
   }

   // random initial guess
   valarray<T> v = rand_vec<T>(A.n,-1.,1.);
   valarray<T> resid = f-A*v;

   // number of iterations to do
   if ( num_itr == -1 ) { // just do normal solve until we reach tol or max_itr
   
      // do iterations
      // norm(...,0) is infinity norm
      num_its_done = 0;
      while ( norm(resid,0) > resid_tol && num_its_done < max_itr ) {
         
         gauss_seidel_it(A,f,v);
   
         resid = f-A*v;
         ++num_its_done;
      }
   }

   else if ( num_itr >= 0 ) { // do exactly num_itr iterations
      // do iterations
      // norm(...,0) is infinity norm
      num_its_todo = static_cast<unsigned>(num_itr);
      num_its_done = 0;
      while ( num_its_done < num_its_todo ) {
         
         gauss_seidel_it(A,f,v);
   
         ++num_its_done;
      }
   }
   
   else {
      cerr << "error: classical_solvers:gauss_seidel: bad input number of "
           << "iterations" << endl;
      exit(-1);
   }

   num_itr = num_its_done;
   return v;
}

template<typename T>
void gauss_seidel_it(const matrix_crs<T>& A, const valarray<T>& f,
                     valarray<T>& v) {
   // Do one iteration of GS
   T LpUv, ajj=0.;

   // sweep down rows of A
   for (size_t row=0; row < A.m; ++row) {
      LpUv = 0.;
   
      // sweep across column.  accumulated (L+U)*v0 and find a_{jj}
      for (size_t ptr=A.row_ptr[row]; ptr < A.row_ptr[row+1]; ++ptr) {
         if ( A.col_ind[ptr] != row ) {// (L+U)*v
            LpUv -= A.val[ptr] * v[A.col_ind[ptr]];
         }
         else {// get a_{jj}
            ajj = A.val[ptr];
         }
      }
   
      // assemble for the update
      v[row] = (LpUv+f[row])/ajj;
   }
}

// }}}

////////////////////////////
// Red-Black Gauss-Seidel //
////////////////////////////
// {{{
template<typename T>
valarray<T> rbgauss_seidel(const matrix_crs<T>& A, const valarray<T>& f,
                         int num_itr) {
   // Red-Black Gauss-Seidel method driver
   // pass in the matrix, RHS vector
   // weight is optional.  defaults to 2./3.
   // num_itr is the number of iterations to do.  defaults to -1, which will
   //   run until we reach the maximum number of iterations or reach tolerance
   //   On return, num_itr will have the number of iterations completed.

   unsigned max_itr = 100;
   T resid_tol = static_cast<T>(10e-8);
   
   unsigned num_its_done, num_its_todo;
   
   // check sizes
   if ( A.m != f.size() ) {
      cerr << "error: classical_solvers:rbgauss_seidel: dimension mismatch"
           << endl;
      exit(-1);
   }

   // it is assumed that A is odd x odd, otherwise RB GS won't work
   if ( A.m % 2 != 1 ) {
      cerr << "error: classical_solvers:rbgauss_seidel: A has even dimensions"
           << endl;
      exit(-1);
   }

   // random initial guess
   valarray<T> v = rand_vec<T>(A.n,-1.,1.);
   valarray<T> resid = f-A*v;

   // number of iterations to do
   if ( num_itr == -1 ) { // just do normal solve until we reach tol or max_itr
   
      // do iterations
      // norm(...,0) is infinity norm
      num_its_done = 0;
      while ( norm(resid,0) > resid_tol && num_its_done < max_itr ) {
         
         rbgauss_seidel_it(A,f,v);
   
         resid = f-A*v;
         ++num_its_done;
      }
   }

   else if ( num_itr >= 0 ) { // do exactly num_itr iterations
      // do iterations
      // norm(...,0) is infinity norm
      num_its_todo = static_cast<unsigned>(num_itr);
      num_its_done = 0;
      while ( num_its_done < num_its_todo ) {
         
         rbgauss_seidel_it(A,f,v);
   
         ++num_its_done;
      }
   }
   
   else {
      cerr << "error: classical_solvers:rbgauss_seidel: bad input number of "
           << "iterations" << endl;
      exit(-1);
   }

   num_itr = num_its_done;
   return v;
}

template<typename T>
void rbgauss_seidel_it(const matrix_crs<T>& A, const valarray<T>& f,
                     valarray<T>& v) {
   // Do one iteration of red-black GS
   T LpUv, ajj=0.;

   // Do the ``red'' updates
   // These are the even numbered nodes for the original problem, which
   // includes the boundary.  The model problem A doesn't include the boundary,
   // so we will actually work on the odd nodes
   for (size_t row=1; row < A.m; row += 2) {
      //cout << "red row = " << row << endl;
      LpUv = 0.;
   
      // sweep across column.  accumulated (L+U)*v0 and find a_{jj}
      for (size_t ptr=A.row_ptr[row]; ptr < A.row_ptr[row+1]; ++ptr) {
         // for red nodes (which uses black nodes)
         if ( A.col_ind[ptr] != row && A.col_ind[ptr] % 2 == 0) {
            //cout << "using col = " << A.col_ind[ptr] << endl;
            LpUv -= A.val[ptr] * v[A.col_ind[ptr]];
         }
         else {// get a_{jj}
            ajj = A.val[ptr];
         }
      }
   
      // assemble for the update
      v[row] = (LpUv+f[row])/ajj;
   }

   // Do the ``black'' updates
   // These are the odd numbered nodes for the original problem, which
   // includes the boundary.  The model problem A doesn't include the boundary,
   // so we will actually work on the even nodes
   for (size_t row=0; row < A.m; row += 2) {
      //cout << "black row = " << row << endl;
      LpUv = 0.;
   
      // sweep across column.  accumulated (L+U)*v0 and find a_{jj}
      for (size_t ptr=A.row_ptr[row]; ptr < A.row_ptr[row+1]; ++ptr) {
         // for black nodes (which uses red nodes)
         if ( A.col_ind[ptr] != row && A.col_ind[ptr] % 2 == 1) {
            LpUv -= A.val[ptr] * v[A.col_ind[ptr]];
         }
         else {// get a_{jj}
            ajj = A.val[ptr];
         }
      }
   
      // assemble for the update
      v[row] = (LpUv+f[row])/ajj;
   }

}

// }}}


// Force instantiation
template valarray<double> wjacobi<double>(const matrix_crs<double>&,
      const valarray<double>&, const double, const int);

template void wjacobi_it<double>(const matrix_crs<double>&,
      const valarray<double>&, const valarray<double>&, valarray<double>&,
      const double);

template valarray<double> gauss_seidel<double>(const matrix_crs<double>&,
      const valarray<double>&, const int);

template void gauss_seidel_it<double>(const matrix_crs<double>&,
      const valarray<double>&, valarray<double>&);

template valarray<double> rbgauss_seidel<double>(const matrix_crs<double>&,
      const valarray<double>&, const int);

template void rbgauss_seidel_it<double>(const matrix_crs<double>&,
      const valarray<double>&, valarray<double>&);

