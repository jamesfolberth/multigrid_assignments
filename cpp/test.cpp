// test.cpp
//
// James Folberth
// Spring 2015

#include <iostream>
#include <iomanip>
#include <valarray>
#include <vector>
#include <random> // randomly generated tests

#include "matrix_coo.hpp"
#include "matrix_crs.hpp"
#include "model_problems.hpp"
#include "classical_solvers.hpp"
#include "multigrid.hpp"

using namespace std;

// Test helpers
template<typename T>
matrix_coo<T> rand_coo_rand_size(void) {
   // {{{
   random_device rd;
   std::mt19937 e2(rd());
   uniform_int_distribution<int> size_dist(1, 9);
   uniform_real_distribution<T> val_dist(-99,99);

   int m = size_dist(e2);
   uniform_int_distribution<int> rind_dist(0,m);
   int n = size_dist(e2);
   uniform_int_distribution<int> cind_dist(0,n);

   uniform_int_distribution<int> nnz_dist(0,(m+1)*(n+1));
   int nnz = nnz_dist(e2);

   vector<unsigned> rind, cind;
   vector<T> vals;

   for (int i=0; i< nnz; ++i) {
      rind.push_back(rind_dist(e2));
      cind.push_back(cind_dist(e2));
      vals.push_back(val_dist(e2));
   }

   return matrix_coo<T>(rind, cind, vals, m, n);
   // }}}
}

template<typename T>
matrix_crs<T> rand_crs_rand_size(void) {
   // {{{
   random_device rd;
   std::mt19937 e2(rd());
   uniform_int_distribution<int> size_dist(1, 9);
   uniform_real_distribution<T> val_dist(-99,99);

   int m = size_dist(e2);
   uniform_int_distribution<int> rind_dist(0,m);
   int n = size_dist(e2);
   uniform_int_distribution<int> cind_dist(0,n);

   uniform_int_distribution<int> nnz_dist(0,(m+1)*(n+1));
   int nnz = nnz_dist(e2);

   vector<unsigned> rind, cind;
   vector<T> vals;

   for (int i=0; i< nnz; ++i) {
      rind.push_back(rind_dist(e2));
      cind.push_back(cind_dist(e2));
      vals.push_back(val_dist(e2));
   }

   return matrix_crs<T>(rind, cind, vals, m, n);
   // }}}
}

template<typename T>
matrix_coo<T> rand_coo(unsigned m, unsigned n, unsigned nnz=0) {
   // {{{
   random_device rd;
   std::mt19937 e2(rd());
   uniform_real_distribution<T> val_dist(-99,99);

   uniform_int_distribution<int> rind_dist(0,m-1);
   uniform_int_distribution<int> cind_dist(0,n-1);

   if ( nnz == 0 ) {
      uniform_int_distribution<int> nnz_dist(0,m*n);
      //uniform_int_distribution<int> nnz_dist(0,m*n/10);
      nnz = nnz_dist(e2);
   }

   vector<unsigned> rind, cind;
   vector<T> vals;

   for (unsigned i=0; i< nnz; ++i) {
      rind.push_back(rind_dist(e2));
      cind.push_back(cind_dist(e2));
      vals.push_back(val_dist(e2));
   }

   return matrix_coo<T>(rind, cind, vals, m, n);
   // }}}
}

template<typename T>
matrix_crs<T> rand_crs(unsigned m, unsigned n, unsigned nnz=0) {
   // {{{
   random_device rd;
   std::mt19937 e2(rd());
   uniform_real_distribution<T> val_dist(-99,99);

   uniform_int_distribution<int> rind_dist(0,m-1);
   uniform_int_distribution<int> cind_dist(0,n-1);

   if ( nnz == 0 ) {
      uniform_int_distribution<int> nnz_dist(0,m*n);
      //uniform_int_distribution<int> nnz_dist(0,m*n/10);
      nnz = nnz_dist(e2);
   }

   vector<unsigned> rind, cind;
   vector<T> vals;

   for (unsigned i=0; i< nnz; ++i) {
      rind.push_back(rind_dist(e2));
      cind.push_back(cind_dist(e2));
      vals.push_back(val_dist(e2));
   }

   return matrix_crs<T>(rind, cind, vals, m, n);
   // }}}
}


///////////////////////
// COO test routines //
///////////////////////
void test_coo_matrix(void) {
   // {{{
   
   // Basic test
   //vector<unsigned> rind = {0,1,2,3};
   //vector<unsigned> cind = {0,1,2,3};
   //vector<double> vals = {1,2,3,4};
   //matrix_coo<double> coomat(rind, cind, vals);
   //cout << coomat << endl;


   // zero rows in beginning, inside, and at end
   //unsigned m=10;
   //unsigned n=5;
   //vector<unsigned> rind = {1,1,1,2,2,4};
   //vector<unsigned> cind = {0,2,3,1,2,3};
   //vector<double> vals = {1,2,3,4,4,5};
   //matrix_coo<double> coomat(rind, cind, vals, m, n);
   //cout << coomat << endl;


   // Another test
   //vector<unsigned> rind, cind;
   //vector<double> vals;
   //unsigned m = 7;
   //rind.resize(m);
   //cind.resize(m);
   //vals.resize(m);

   //for (unsigned i=0; i<m; ++i) {
   //   rind[i] = i;
   //   cind[i] = i;
   //   vals[i] = 1e0*sqrt(double(i));
   //}
   //matrix_coo<double> coomat(rind, cind, vals);
   //cout << coomat << endl;


   // Build random mat and print
   //matrix_coo<double> A = rand_coo_rand_size<double>();
   matrix_coo<double> A = rand_coo<double>(4,5);
   cout << A << endl;
   A.print_full();
   
   // }}}
}

void test_coo_eye(void) {
   // {{{

   cout << "eye_coo(7,7) = " << endl;
   matrix_coo<double> I = eye_coo<double>(7,7);
   I.print_full(); 

   cout << "eye_coo(5,7) = " << endl;
   I = eye_coo<double>(5,7);
   I.print_full();

   cout << "eye_coo(7,5) = " << endl;
   I = eye_coo<double>(7,5);
   I.print_full();

   // }}}
}

void test_coo_scalar(void) {
   // {{{
   matrix_coo<double> A = rand_coo<double>(8,5);
   matrix_coo<double> B = A;

   cout << "A = " << endl;
   B.print_full();

   cout << "A *= 2.5" << endl;
   B *= 2.5; B.print_full(); B = A;

   cout << "2.5*A" << endl;
   (2.5*A).print_full();

   cout << "A*2.5" << endl;
   (A*2.5).print_full();

   cout << "A *= 0.0 then A.clean()" << endl;
   B *= 0; B.clean(); cout << B << endl; B = A;

   cout << "A /= 2.5" << endl;
   B /= 2.5; B.print_full(); B = A;

   cout << "A/2.5" << endl;
   (A/2.5).print_full();
 
   // }}}
}

void test_coo_add(void) {
   // {{{
   matrix_coo<double> A = rand_coo<double>(3,5);
   matrix_coo<double> B = rand_coo<double>(3,5);
   //matrix_crs<double> B = A; B *= -1.0;

   cout << "A = " << endl;
   //cout << A << endl;
   A.print_full();

   cout << "B = " << endl;
   //cout << B << endl;
   B.print_full();

   cout << "C = A + B" << endl;
   matrix_coo<double> C = A + B;
   //cout << A << endl;
   C.print_full();

   cout << "A += B" << endl;
   C = A; C += B;
   //cout << C << endl;
   C.print_full();

   cout << "A -= B" << endl;
   C = A; C -= B;
   //cout << C << endl;
   C.print_full();

   // }}}
}



///////////////////////
// CRS test routines //
///////////////////////
void test_crs_matrix(void) {
   // {{{

   // Basic test
   //vector<unsigned> rind = {0,1,2,3};
   //vector<unsigned> cind = {0,1,2,3};
   //vector<double> vals = {1,2,3,4};
   //matrix_crs<double> crsmat(rind, cind, vals);
   //cout << crsmat << endl;


   // zero rows in beginning, inside, and at end
   //unsigned m=10;
   //unsigned n=5;
   //vector<unsigned> rind = {1,1,1,2,2,4};
   //vector<unsigned> cind = {0,2,3,1,2,3};
   //vector<double> vals = {1,2,3,4,4,5};
   //matrix_crs<double> crsmat(rind, cind, vals, m, n);
   //cout << crsmat << endl;


   // Another test
   //vector<unsigned> rind, cind;
   //vector<double> vals;
   //unsigned m = 4;
   //rind.resize(m);
   //cind.resize(m);
   //vals.resize(m);

   //for (unsigned i=0; i<m; ++i) {
   //   rind[i] = i;
   //   cind[i] = i;
   //   vals[i] = 1e0*sqrt(double(i+1));
   //}
   //matrix_crs<double> crsmat(rind, cind, vals);
   //cout << crsmat << endl;
   //crsmat.print_full();

   // Build random mat
   //matrix_crs<double> A = rand_crs_rand_size<double>();
   matrix_crs<double> A = rand_crs<double>(4,5);
   cout << A << endl;
   A.print_full();

   // }}}
}

void test_crs_eye(void) {
   // {{{

   cout << "eye_crs(7,7) = " << endl;
   matrix_crs<double> I = eye_crs<double>(7,7);
   I.print_full(); 

   cout << "eye_crs(5,7) = " << endl;
   I = eye_crs<double>(5,7);
   I.print_full();

   cout << "eye_crs(7,5) = " << endl;
   I = eye_crs<double>(7,5);
   I.print_full();

   // }}}
}

void test_crs_scalar(void) {
   // {{{
   matrix_crs<double> A = rand_crs<double>(8,5);
   matrix_crs<double> B = A;

   cout << "A = " << endl;
   B.print_full();

   cout << "A *= 2.5" << endl;
   B *= 2.5; B.print_full(); B = A;

   cout << "2.5*A" << endl;
   (2.5*A).print_full();

   cout << "A*2.5" << endl;
   (A*2.5).print_full();

   cout << "A *= 0.0 then A.clean()" << endl;
   B *= 0; B.clean(); cout << B << endl; B = A;

   cout << "A /= 2.5" << endl;
   B /= 2.5; B.print_full(); B = A;

   cout << "A/2.5" << endl;
   (A/2.5).print_full();

   // }}}
}

void test_crs_add(void) {
   // {{{
   matrix_crs<double> A = rand_crs<double>(3,5);
   matrix_crs<double> B = rand_crs<double>(3,5);
   //matrix_crs<double> B = A; B *= -1.0;

   cout << "A = " << endl;
   //cout << A << endl;
   A.print_full();

   cout << "B = " << endl;
   //cout << B << endl;
   B.print_full();

   cout << "C = A + B" << endl;
   matrix_crs<double> C = A + B;
   //cout << C << endl;
   C.print_full();

   cout << "A += B" << endl;
   C = A; C += B;
   //cout << C << endl;
   C.print_full();

   cout << "A -= B" << endl;
   C = A; C -= B;
   //cout << C << endl;
   C.print_full();

   // }}}
}

void test_crs_kron(void) {
   // {{{
   matrix_crs<double> A = rand_crs<double>(2,2);
   matrix_crs<double> B = rand_crs<double>(2,3);
   matrix_crs<double> C = kron<double>(A,B);

   cout << "A = " << endl;
   A.print_full();

   cout << "B = " << endl;
   B.print_full();

   cout << "kron(A,B) = " << endl;
   //cout << C << endl;
   C.print_full();
   // }}}
}

void test_crs_matvec(void) {
   // {{{
   unsigned m = 5, n=6;
   matrix_crs<double> A = rand_crs<double>(m,n);
   valarray<double> v(1.,n),u;

   cout << "A = " << endl;
   A.print_full();

   cout << "v = " << endl;
   print_vector(v);

   cout << "u = " << endl;
   u = A*v;
   print_vector(u);
   
   // }}}
}


////////////////////
// Model problems //
////////////////////
void test_model_problems(void) {
   // {{{
   //matrix_crs<double> A1 = model_problem_1d<double>(3,1.);
   //cout << A1 << endl;
   //A1.print_full();

   matrix_crs<double> A2 = model_problem_2d<double>(2,3,1.);
   cout << A2 << endl;
   A2.print_full();
   // }}}
}


///////////////////////
// Classical Solvers //
///////////////////////
void test_wjacobi(void) {
   // {{{
   unsigned Lx = 3;
   matrix_crs<double> A = model_problem_1d(Lx,0.);
   valarray<double> f(0.,pow(2,Lx)-1);
   valarray<double> v0,v1;
   int zero = 0;

   // use the driver to set things up
   v0 = wjacobi<double>(A,f,2./3.,zero);
   v0 /= norm(v0,0);
   v1 = v0;

   for (unsigned i=0; i < 100; ++i) {
      v0 = v1;
      //wjacobi_it<double>(A,f,v0,v1,2./3.);
      wjacobi_ip<double>(A,f,v1,1,2./3.);

      cout << "\\|error\\|_inf = " << norm(v1,0) << endl;
   }
   
   // }}}
}

void test_gauss_seidel(void) {
   // {{{
   unsigned Lx = 6;
   matrix_crs<double> A = model_problem_1d(Lx,0.);
   valarray<double> f(0.,pow(2,Lx)-1);
   valarray<double> v;
   int zero = 0;

   // use the driver to set things up
   v = gauss_seidel<double>(A,f,zero);
   v /= norm(v,0);

   for (unsigned i=0; i < 100; ++i) {
      gauss_seidel_it<double>(A,f,v);

      cout << "\\|error\\|_inf = " << norm(v,0) << endl;
   }
   
   // }}}
}

void test_rbgauss_seidel(void) {
   // {{{
   unsigned Lx = 6;
   matrix_crs<double> A = model_problem_1d(Lx,0.);
   valarray<double> f(0.,pow(2,Lx)-1);
   valarray<double> v;
   int zero = 0;

   // use the driver to set things up
   v = rbgauss_seidel<double>(A,f,zero);
   v /= norm(v,0);

   for (unsigned i=0; i < 100; ++i) {
      //rbgauss_seidel_it<double>(A,f,v);
      rbgauss_seidel_ip<double>(A,f,v,1);

      cout << "\\|error\\|_inf = " << norm(v,0) << endl;
   }
   
   // }}}
}


///////////////
// Multigrid //
///////////////
void test_mg_1d_intergrid_operators(void) {
   // {{{
   int L = 3;

   matrix_crs<double> P = operator_1d_interp_lin<double>(L);
   cout << "P = " << endl;
   P.print_full();
  
   matrix_crs<double> Rinj = operator_1d_restrict_inj<double>(L);
   cout << "Rinj = " << endl;
   Rinj.print_full();

   matrix_crs<double> Rfull = operator_1d_restrict_full<double>(L);
   cout << "Rfull = " << endl;
   Rfull.print_full();
   // }}}
}

void test_mg_1d_vcycle(void) {
   // {{{ 
   unsigned L = 17, Lx = 20, n = pow(2,Lx)-1;
   double sigma = 0., resid_nrm, err_nrm, resid_nrm_prev;
   valarray<double> f(0.,n), u(0.,n), v0(0.,n), v(0.,n),
      resid(0.,n), err(0.,n);

   double xi = 0.,C = 2., k = 3.;

   int mode = 2;
   switch (mode) {
      case 0:
         // zero RHS
         v0 = rand_vec<double>(n,-1.,1.);
         f *= 0.; u *= 0.;
         break;

      case 1:
         // sin(pi*k*x) RHS
         v0 = rand_vec<double>(n,-1.,1.);

         for (unsigned i=0; i < n; ++i) {
            f[i] = C*sin(k*_PI_*double(i+1-0)/double(n+1));
            u[i] = C/(pow(_PI_*k,2.)+sigma)*sin(k*_PI_*double(i+1-0)/double(n+1));
            //cout << setprecision(10) << double(i+1-0)/double(n+1) << endl;
         }

         // need to account for dirichlet BCs
         f[0] += pow(n+1,2)*0.;
         f[n-1] += pow(n+1,2)*0.;

         break;

      case 2:
         // sin(pi*x)*exp(-x) RHS
         v0 = rand_vec<double>(n,-1.,1.);

         for (unsigned i=0; i < n; ++i) {
            xi = double(i+1-0)/double(n+1);
            u[i] = sin(_PI_*xi)*exp(-xi);
            f[i] = exp(-xi)*((-1+_PI_*_PI_+sigma)*sin(_PI_*xi) 
                  + 2.*_PI_*cos(_PI_*xi));
         }

         // need to account for dirichlet BCs
         f[0] += pow(n+1,2)*0.;
         f[n-1] += pow(n+1,2)*0.;

         break;
 
      default:
         cerr << "error: test.cpp:test_mg_1d_vcycle: bad mode" << endl;
         exit(-1);
   }

  
   // make a function to build A that depends only on the grid level
   // (so fix sigma in the beginning)
   auto build_A = [=](unsigned _Lx) -> matrix_crs<double>  // [=] pass by copy
      {return model_problem_1d<double>(_Lx,sigma);};
   auto build_P = [](unsigned _Lx) -> matrix_crs<double>
      {return operator_1d_interp_lin<double>(_Lx);};
   auto build_R = [](unsigned _Lx) -> matrix_crs<double>
      //{return operator_1d_restrict_inj<double>(_Lx);};
      {return operator_1d_restrict_full<double>(_Lx);};

   auto smoother = [](const matrix_crs<double>& _A, const valarray<double>& _f,
         valarray<double>& _v, unsigned _num_itr)
      {wjacobi_ip<double>(_A,_f,_v,_num_itr);};
      //{gauss_seidel_ip<double>(_A,_f,_v,_num_itr);};
      //{rbgauss_seidel_ip<double>(_A,_f,_v,_num_itr);};


   // build levels
   vector<level<double>> levels = build_levels<double>(build_A, f, build_P,
         build_R, smoother, L, Lx, v0);

   // do V cycle
   cout << "Matrix A is " << levels[0].A.m << " x " << levels[0].A.n << endl;
   
   resid = levels[0].f - levels[0].A*levels[0].v;
   err = u-levels[0].v;
   resid_nrm_prev = pow(levels[0].A.n+1,-0.5)*norm(resid,2);
   err_nrm = pow(levels[0].A.n+1,-0.5)*norm(err,2);
   
   cout << "||r||_h = " << _PRINT_VECTOR_FORMAT_ << resid_nrm_prev
        << "  ||e||_h = " << _PRINT_VECTOR_FORMAT_ << err_nrm << endl;
 
   for (unsigned i = 0; i < 10; ++i) {
      vcycle(levels, levels.begin(), 2, 1);

      resid = levels[0].f - levels[0].A*levels[0].v;
      err = u-levels[0].v;
      resid_nrm = pow(levels[0].A.n+1,-0.5)*norm(resid,2);
      err_nrm = pow(levels[0].A.n+1,-0.5)*norm(err,2);

      cout << "||r||_h = " << _PRINT_VECTOR_FORMAT_ << resid_nrm
           << "  ||e||_h = " << _PRINT_VECTOR_FORMAT_ << err_nrm
           << "  ratio = " << _PRINT_VECTOR_FORMAT_ 
           << resid_nrm/resid_nrm_prev << "  at i = " << i << endl;

      resid_nrm_prev = resid_nrm;
   }

   // }}}
}



int main() {
  
   // COO
   //test_coo_matrix();
   //test_coo_eye();
   //test_coo_scalar();
   //test_coo_add();

   // CRS
   //test_crs_matrix();
   //test_crs_eye();
   //test_crs_scalar();
   //test_crs_add();
   //test_crs_kron();
   //test_crs_matvec();
   
   // Model problems
   //test_model_problems();

   // Classical solvers
   //test_wjacobi();
   //test_gauss_seidel();
   //test_rbgauss_seidel();
   
   // Multigrid
   //test_mg_1d_intergrid_operators();
   test_mg_1d_vcycle();

   return 0;
}
