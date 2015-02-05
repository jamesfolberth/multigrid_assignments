// utils.cpp
//
// James Folberth
// Spring 2015

#include <iostream>

#include "utils.hpp"

using namespace std;


//////////
// Misc //
//////////

// Integer pow
int pow(int b, int e) {
   if ( e == 0 ) return 1;
   else if ( e == 1 ) return b;
   else if ( e < 0 ) {
      cerr << "utils.cpp:pow: negative exponent not handled" << endl;
      exit(-1);
   }
   else return pow(b, e-1);
}




