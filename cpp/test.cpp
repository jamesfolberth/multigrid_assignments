// test.cpp
//
// James Folberth
// Spring 2015

#include <iostream>
#include <vector>

#include "matrix_coo.hpp"

using namespace std;

void test_matrix_coo(void) {

   vector<unsigned> rind, cind;
   vector<double> vals;

   unsigned m = 10;
   rind.resize(m);
   cind.resize(m);
   vals.resize(m);

   for (unsigned i=0; i<m; ++i) {
      rind[i] = i/2;
      cind[i] = 10-i;
      vals[i] = double(i);
   }

   matrix_coo<double> coomat(rind, cind, vals, m, m);
   cout << coomat << endl;

}

void test_matrix_crs(void) {

   vector<unsigned> rind, cind;
   vector<double> vals;

   unsigned m = 10;
   rind.resize(m);
   cind.resize(m);
   vals.resize(m);

   for (unsigned i=0; i<m; ++i) {
      rind[i] = i;
      cind[i] = i;
      vals[i] = double(i);
   }

   //matrix_crs<double> crsmat(rind, cind, vals, m, m);
   //cout << crsmat << endl;

}

int main() {

   test_matrix_coo();

   return 0;
}
