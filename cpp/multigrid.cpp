// multigrid.cpp
//
// James Folberth
// Spring 2015

#include "multigrid.hpp"

using namespace std;

///////////
// Level //
///////////
// {{{
template<typename T>
vector<level<T>> build_levels(function<matrix_crs<T>(unsigned)> build_A,
      valarray<T> f, function<matrix_crs<T>(unsigned)> build_P,
      function<matrix_crs<T>(unsigned)> build_R,
      function<void(const matrix_crs<T>&, const valarray<T>&,
         valarray<T>&, unsigned)> smoother_ip,
      unsigned L, unsigned Lx, const valarray<T>& v0) {
   // Set up A,P,R,smoother, etc. for levels.  This is just a generic setup
   // for multigrid
   //
   // build_A, build_P, and build_R take exactly one argument: Lx
   // Use a lambda function to define build_A if you're using something
   // like model_problem_1d(Lx, sigma)

   // build and populate vector of levels
   vector<level<T>> levels(L);

   // A,P,R
   for (unsigned l = 0; l < L; ++l) {
      levels[l].A = build_A(Lx-l);
      levels[l].f = valarray<T>(0.,levels[l].A.m);
      levels[l].v = valarray<T>(0.,levels[l].A.n);
      levels[l].P = build_P(Lx-l);
      levels[l].R = build_R(Lx-l);
      levels[l].smoother_ip = smoother_ip;
   }

   // initial RHS and v
   levels[0].f = f;
   levels[0].v = v0;

   //for (auto it = levels.begin(); it != levels.end(); ++it) {
   //   //it->A.print_full();
   //   //it->P.print_full();
   //   it->R.print_full();
   //}

   return levels;
}

template<typename T>
vector<level<T>> build_levels(function<matrix_crs<T>(unsigned)> build_A,
      valarray<T> f, function<matrix_crs<T>(unsigned)> build_P,
      function<matrix_crs<T>(unsigned)> build_R,
      function<void(const matrix_crs<T>&, const valarray<T>&,
         valarray<T>&, unsigned)> smoother_ip,
      unsigned L, unsigned Lx) {
   // Set up A,P,R,smoother, etc. for levels.  This is just a generic setup
   // for multigrid.  
   //
   // Generate a random vector v0, as opposed to the caller supplying one

   return build_levels(build_A, f, build_P, build_R, smoother_ip, L, Lx,
         rand_vec<T>(pow(2,Lx)-1, static_cast<T>(-1.), static_cast<T>(1.)));
}
// }}}


////////////
// Cycles //
////////////
// {{{
template<typename T>
void vcycle(vector<level<T>>& levels, typename vector<level<T>>::iterator it,
      unsigned nu1, unsigned nu2) {                                      
   // Perform a (nu1,nu2) V cycle on the levels                        
   // the iterator it points to our current position in the levels vector.

   // pre-smooth
   it->smoother_ip(it->A, it->f, it->v, nu1);
   cout << "pre-smooth error (" << it->v.size() << ") = " << norm(it->v,0) << endl;

   if ( it == levels.end()-1 ) { // if we're on the coarsest grid
      // smooth the dick out of it
      // TODO direct solve
      cout << "Direct solve start (" << it->v.size() << ") = " << norm(it->v,0) << endl;
      it->smoother_ip(it->A, it->f, it->v, 100);
      cout << "Direct solve error (" << it->v.size() << ") = " << norm(it->v,0) << endl;
   }

   else { // we're not on the coarsest grid

      // prepare to coarsen
      valarray<T> temp = it->f;
      temp -= (it->A)*(it->v);

      //print_vector((it->f) - (it->A)*(it->v)); // TODO why doesn't this work?
      next(it)->f = (it->R)*temp;
      next(it)->v *= 0.;

      // move to coarser grid
      vcycle(levels, next(it), nu1, nu2);

      // correct this grid using coarse grid
      it->v += (next(it)->P)*(next(it)->v);
   }

   // post-smooth
   it->smoother_ip(it->A, it->f, it->v, nu2);
   cout << "post-smooth error (" << it->v.size() << ") = " << norm(it->v,0) << endl;
   return;
}


// }}}


/////////////////////////
// Intergrid Operators //
/////////////////////////
// {{{
template<typename T>
matrix_crs<T> operator_1d_interp_lin(unsigned l) {
// form the linear operator P that performs linear interpolation from grid
// l to grid l+1.  Number of points on the initial grid is 2^l-1

   unsigned m = pow(2,l+1)-1, n = pow(2,l)-1;

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
template vector<level<double>> build_levels(
      function<matrix_crs<double>(unsigned)> build_A, valarray<double> f,
      function<matrix_crs<double>(unsigned)> build_P,
      function<matrix_crs<double>(unsigned)> build_R,
      function<void(const matrix_crs<double>&, 
         const valarray<double>&, valarray<double>&, unsigned)>,
      unsigned L, unsigned Lx, const valarray<double>& v0);

template vector<level<double>> build_levels(
      function<matrix_crs<double>(unsigned)> build_A, valarray<double> f,
      function<matrix_crs<double>(unsigned)> build_P,
      function<matrix_crs<double>(unsigned)> build_R,
      function<void(const matrix_crs<double>&, 
         const valarray<double>&, valarray<double>&, unsigned)>,
      unsigned L, unsigned Lx);


template void vcycle(vector<level<double>>& levels, 
      vector<level<double>>::iterator, unsigned nu1, unsigned nu2);

template matrix_crs<double> operator_1d_interp_lin<double>(unsigned);
template matrix_crs<double> operator_1d_restrict_inj<double>(unsigned);
template matrix_crs<double> operator_1d_restrict_full<double>(unsigned);


