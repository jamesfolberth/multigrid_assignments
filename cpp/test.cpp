// test.cpp
//
// James Folberth
// Spring 2015

#include <iostream>
#include <vector>
#include <random> // randomly generated tests

#include "matrix_coo.hpp"
#include "matrix_crs.hpp"

using namespace std;

void test_matrix_coo(void) {
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
   //unsigned m = 10;
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


   // Build random mat
   random_device rd;
   std::mt19937 e2(rd());
   uniform_int_distribution<int> size_dist(1, 9);
   uniform_real_distribution<double> val_dist(-99,99);


   int m = size_dist(e2);
   uniform_int_distribution<int> rind_dist(0,m);
   int n = size_dist(e2);
   uniform_int_distribution<int> cind_dist(0,n);

   uniform_int_distribution<int> nnz_dist(0,(m+1)*(n+1));
   int nnz = nnz_dist(e2);

   vector<unsigned> rind, cind;
   vector<double> vals;

   for (int i=0; i< nnz; ++i) {
      rind.push_back(rind_dist(e2));
      cind.push_back(cind_dist(e2));
      vals.push_back(val_dist(e2));
   }

   matrix_coo<double> coomat(rind, cind, vals);
   cout << coomat << endl;

   coomat.print_full();

   // }}}
}

void test_matrix_crs(void) {
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
   //unsigned m = 10;
   //rind.resize(m);
   //cind.resize(m);
   //vals.resize(m);

   //for (unsigned i=0; i<m; ++i) {
   //   rind[i] = i;
   //   cind[i] = i;
   //   vals[i] = 1e0*sqrt(double(i));
   //}
   //matrix_crs<double> crsmat(rind, cind, vals);
   //cout << crsmat << endl;

   // Build random mat
   random_device rd;
   std::mt19937 e2(rd());
   uniform_int_distribution<int> size_dist(1, 9);
   uniform_real_distribution<double> val_dist(-99,99);


   int m = size_dist(e2);
   uniform_int_distribution<int> rind_dist(0,m);
   int n = size_dist(e2);
   uniform_int_distribution<int> cind_dist(0,n);

   uniform_int_distribution<int> nnz_dist(0,(m+1)*(n+1));
   int nnz = nnz_dist(e2);

   vector<unsigned> rind, cind;
   vector<double> vals;

   for (int i=0; i< nnz; ++i) {
      rind.push_back(rind_dist(e2));
      cind.push_back(cind_dist(e2));
      vals.push_back(val_dist(e2));
   }

   matrix_crs<double> crsmat(rind, cind, vals);
   cout << crsmat << endl;

   crsmat.print_full();

   // }}}
}

int main() {
   
   //test_matrix_coo();
   test_matrix_crs();

   return 0;
}
