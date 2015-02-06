// test.cpp
//
// James Folberth
// Spring 2015

#include <iostream>
#include <vector>
#include <random> // randomly generated tests

#include "matrix_coo.hpp"
#include "matrix_crs.hpp"
#include "model_problems.hpp"

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
matrix_coo<T> rand_coo(unsigned m, unsigned n) {
   // {{{
   random_device rd;
   std::mt19937 e2(rd());
   uniform_real_distribution<T> val_dist(-99,99);


   uniform_int_distribution<int> rind_dist(0,m-1);
   uniform_int_distribution<int> cind_dist(0,n-1);

   uniform_int_distribution<int> nnz_dist(0,m*n);
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
matrix_crs<T> rand_crs(unsigned m, unsigned n) {
   // {{{
   random_device rd;
   std::mt19937 e2(rd());
   uniform_real_distribution<T> val_dist(-99,99);


   uniform_int_distribution<int> rind_dist(0,m-1);
   uniform_int_distribution<int> cind_dist(0,n-1);

   uniform_int_distribution<int> nnz_dist(0,m*n);
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

   cout << "A += 2.5" << endl;
   B += 2.5; B.print_full(); B = A;

   cout << "A -= 2.5" << endl;
   B -= 2.5; B.print_full(); B = A;

   cout << "A *= 2.5" << endl;
   B *= 2.5; B.print_full(); B = A;

   cout << "A /= 2.5" << endl;
   B /= 2.5; B.print_full(); B = A;

   cout << "A *= 0.0 then A.clean()" << endl;
   B *= 0; B.clean(); B.print_full(); B = A;
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


////////////////////
// Model problems //
////////////////////
void test_model_problems(void) {

   matrix_crs<double> A1 = model_problem_1d<double>(3,1.);
   cout << A1 << endl;
   A1.print_full();

}

int main() {
  
   // COO
   //test_coo_matrix();
   //test_coo_eye();
   //test_coo_scalar();

   // CRS
   //test_crs_matrix();
   //test_crs_eye();
   //test_crs_scalar();
   test_crs_add();
   
   // Model problems
   //test_model_problems();

   return 0;
}
