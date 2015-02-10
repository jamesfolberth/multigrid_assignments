// utils.cpp
//
// James Folberth
// Spring 2015

#include "utils.hpp"

using namespace std;

//////////////////
// Vector stuff //
//////////////////

// print a vector as a column vector
template<typename T>
void print_vector(vector<T>& v) {
   for (auto it = v.begin(); it != v.end(); ++it) {
      cout << "    " << setw(_PRINT_VECTOR_WIDTH_) << setfill(' ') 
           << setprecision(_PRINT_VECTOR_PREC_)
           << static_cast<double>(*it) << endl;
   }
   cout << endl;
}

// return a vector with uniform[0,1] random entries
template<typename T>
vector<T> rand_vec(const unsigned m, const T low, const T high) {
   random_device rd;
   std::mt19937 e2(rd());
   uniform_real_distribution<T> val_dist(low,high);

   vector<T> v;
   v.reserve(m);

   for (unsigned i=0; i<m; ++i) {
      v.push_back(val_dist(e2));
   }
   
   return v;
}

// vector p-norms
// use p=0 as infinity norm
// p=2 is like BLAS dnrm2
template<typename T>
T norm(const vector<T>& v, const unsigned p) {
   T res = 0, scale = 0., absit, ssq = 1., tmp;

   if ( v.size() < 1 ) {
      return static_cast<T>(0.);
   }

   else if ( v.size() == 1 ) {
      return abs(v[0]);
   }

   else {
      switch (p) {
         case 0:
            for (auto it = v.begin(); it != v.end(); ++it) {
               if ( res < abs(*it) )
                  res = abs(*it);
            }
            return res;

         case 1:
            for (auto it = v.begin(); it != v.end(); ++it) {
               res += abs(*it);
            }
            return res;

         case 2:
            for (auto it = v.begin(); it != v.end(); ++it) {
               if ( *it != 0. ) {
                  absit = abs(*it);
                  if ( scale < absit ) {
                     tmp = scale/absit;
                     ssq = 1. + ssq*tmp*tmp;
                     scale = absit;
                  }
                  else {
                     tmp = absit/scale;
                     ssq += tmp*tmp;
                  }
               }
            }
            res = scale*sqrt(ssq);
            return res;

         default:
            cerr << "error: utils:norm(vector): unsupported p-norm: p = " << p 
                 << endl;
            exit(-1);
      }  
   }
}

// operations on vectors
template<typename T>
vector<T> operator*(const vector<T>& a, const T val) {
   vector<T> res;
   res.reserve(a.size());
   for (auto it = a.begin(); it != a.end(); ++it) {
      res.push_back((*it)*val);
   }
   return res;
}


template<typename T>
vector<T> operator+(const vector<T>& a, const vector<T>& b) {
   vector<T> res;

   if ( a.size() != b.size() ) {
      cerr << "error: utils:vector:+: dimension mismatch" << endl;
      exit(-1);
   }

   res.reserve(a.size());
   for (unsigned i=0; i < a.size(); ++i) {
      res.push_back(a[i]+b[i]);
   }

   return res;
}


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

// Force instantiation
template void print_vector<double>(vector<double>& v);
template vector<double> rand_vec<double>(const unsigned,
                                         const double, const double);
template double norm<double>(const vector<double>&, const unsigned);

template vector<double> operator*(const vector<double>&, const double);

template vector<double> operator+(const vector<double>&,
                                  const vector<double>&);
//template vector<double> operator-(const vector<double>&,
//                                  const vector<double>&);



