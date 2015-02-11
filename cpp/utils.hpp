// utils.hpp
//
// James Folberth
// Spring 2015

#ifndef _UTILS_HPP_
#define _UTILS_HPP_

#include <iostream>
#include <iomanip>
#include <random>
#include <valarray>
#include <vector>

using namespace std;

// Debug info
// 0 - no debug info
// 1 - some debug info
// 2 - lots of debug info
#define _DEBUG_ 1


// Output formatting
#define _PRINT_SPARSE_PREC_ 5

#define _PRINT_FULL_PREC_ 4
#define _PRINT_FULL_WIDTH_ 10

#define _PRINT_VECTOR_PREC_ 5
#define _PRINT_VECTOR_WIDTH_ 12


#define _ELEMENT_ZERO_TOL_ 10e-15

// misc
#define MIN(a,b) (((a)<(b)) ? a : b)
#define MAX(a,b) (((a)<(b)) ? b : a)

#define _PI_ 3.14159265358979323846264338327950288419716

//////////////////
// Vector stuff //
//////////////////
template<typename T>
void print_vector(valarray<T>& v);

template<typename T>
valarray<T> rand_vec(const unsigned m, const T low=0., const T high=1.);

template<typename T>
T norm(const valarray<T>& v, const unsigned p);


//// operations on vectors
//template<typename T>
//vector<T> operator*(const vector<T>& a, const T val);
//
//template<typename T>
//vector<T> operator+(const vector<T>& a, const vector<T>& b);
//
//template<typename T>
//vector<T> operator-(const vector<T>& a, const vector<T>& b);




//////////
// Misc //
//////////
int pow(int b, int e);

#endif
