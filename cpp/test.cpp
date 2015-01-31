// test.cpp
//
// James Folberth
// Spring 2015

#include <iostream>
#include <vector>

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
   vector<unsigned> rind, cind;
   vector<double> vals;
   unsigned m = 10;
   rind.resize(m);
   cind.resize(m);
   vals.resize(m);

   for (unsigned i=0; i<m; ++i) {
      rind[i] = i;
      cind[i] = i;
      vals[i] = 1e0*sqrt(double(i));
   }
   matrix_coo<double> coomat(rind, cind, vals);
   cout << coomat << endl;

   matrix_coo<double>* newmat = new matrix_coo<double>(rind, cind, vals);

   //coomat.to_crs();
   //coomat.print_full();

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
   vector<unsigned> rind, cind;
   vector<double> vals;
   unsigned m = 10;
   rind.resize(m);
   cind.resize(m);
   vals.resize(m);

   for (unsigned i=0; i<m; ++i) {
      rind[i] = i;
      cind[i] = i;
      vals[i] = 1e0*sqrt(double(i));
   }
   matrix_crs<double> crsmat(rind, cind, vals);
   cout << crsmat << endl;

   crsmat.print_full();

   // }}}
}

int main() {

   test_matrix_coo();
   //test_matrix_crs();

   return 0;
}
