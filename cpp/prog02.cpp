// prog02.cpp
//
// James Folberth
// Spring 2015

#include <cmath>
#include <iostream>
#include <fstream>
#include <valarray>

#include "matrix_crs.hpp"
#include "model_problems.hpp"
#include "utils.hpp"
#include "classical_solvers.hpp"

using namespace std;

template<typename T>
void wjacobi_write_errors(const matrix_crs<T>& A, const valarray<T>& f,
                    valarray<T>& v0, const T w, const unsigned num_its_todo,
                    ofstream& file) {
   // driver to run iterative method and plot iteration number and error
   // to file for printing with gnuplot

   // print initial error (e = u-v = 0-v = -v) to file 
   valarray<T> v1 = v0;
   file << 0 << "  " << norm(v0,0) << endl;

   // do iterations of method and print errors to file
   for (unsigned i=1; i <= num_its_todo; ++i) {
      v0 = v1;
      wjacobi_it<double>(A,f,v0,v1,2./3.);

      file << i << "  " << norm(v1,0) << endl;
   }
}

template<typename T>
void gauss_seidel_write_errors(const matrix_crs<T>& A, const valarray<T>& f,
                    valarray<T>& v, const unsigned num_its_todo,
                    ofstream& file) {
   // driver to run iterative method and plot iteration number and error
   // to file for printing with gnuplot

   // print initial error (e = u-v = 0-v = -v) to file 
   file << 0 << "  " << norm(v,0) << endl;

   // do iterations of method and print errors to file
   for (unsigned i=1; i <= num_its_todo; ++i) {
      gauss_seidel_it<double>(A,f,v);

      file << i << "  " << norm(v,0) << endl;
   }
}

template<typename T>
void rbgauss_seidel_write_errors(const matrix_crs<T>& A, const valarray<T>& f,
                    valarray<T>& v, const unsigned num_its_todo,
                    ofstream& file) {
   // driver to run iterative method and plot iteration number and error
   // to file for printing with gnuplot

   // print initial error (e = u-v = 0-v = -v) to file 
   file << 0 << "  " << norm(v,0) << endl;

   // do iterations of method and print errors to file
   for (unsigned i=1; i <= num_its_todo; ++i) {
      rbgauss_seidel_it<double>(A,f,v);

      file << i << "  " << norm(v,0) << endl;
   }
}


void make_2_3_plot_data(unsigned method) {
   // driver to call iterative method error writing functions
   // method = 1 -> wjacobi
   // method = 2 -> gauss-seidel
   // method = 3 -> red-black gauss-seidel

   // set up the problem and evecs
   matrix_crs<double> A = model_problem_1d<double>(6,0.0);
   valarray<double> f(0.,A.m); // zero RHS
   //valarray<double> v; // dummy return variable
   valarray<double> v1(0.,A.n), v3(0.,A.n), v6(0.,A.n);

   for (size_t j=0; j < v1.size(); ++j) {
      v1[j] = sin(static_cast<double>(1*(j+1))*_PI_ 
                  / static_cast<double>(A.n+1));
      v3[j] = sin(static_cast<double>(3*(j+1))*_PI_ 
                  / static_cast<double>(A.n+1));
      v6[j] = sin(static_cast<double>(6*(j+1))*_PI_ 
                  / static_cast<double>(A.n+1));
   }

   ofstream k1file, k3file, k6file;

   // set up output streams
   if ( method == 1 ) {
      k1file.open("figures/prog02/2_3_wj_k1_error.txt", ios::out | ios::trunc);
      k3file.open("figures/prog02/2_3_wj_k3_error.txt", ios::out | ios::trunc);
      k6file.open("figures/prog02/2_3_wj_k6_error.txt", ios::out | ios::trunc);

      // run methods and write to files
      if ( k1file.is_open() && k3file.is_open() && k6file.is_open() ) {

         wjacobi_write_errors(A,f,v1,2./3.,100, k1file);
         wjacobi_write_errors(A,f,v3,2./3.,100, k3file);
         wjacobi_write_errors(A,f,v6,2./3.,100, k6file);

         k1file.close();
         k3file.close();
         k6file.close();
      }
      
      else {
         cerr << "file output error" << endl;
         exit(-1);
      }
   }

   else if ( method == 2) {
      k1file.open("figures/prog02/2_3_gs_k1_error.txt", ios::out | ios::trunc);
      k3file.open("figures/prog02/2_3_gs_k3_error.txt", ios::out | ios::trunc);
      k6file.open("figures/prog02/2_3_gs_k6_error.txt", ios::out | ios::trunc);

      // run methods and write to files
      if ( k1file.is_open() && k3file.is_open() && k6file.is_open() ) {

         gauss_seidel_write_errors(A,f,v1,100, k1file);
         gauss_seidel_write_errors(A,f,v3,100, k3file);
         gauss_seidel_write_errors(A,f,v6,100, k6file);

         k1file.close();
         k3file.close();
         k6file.close();
      }
      
      else {
         cerr << "file output error" << endl;
         exit(-1);
      }
   }

   else if ( method == 3) {
      k1file.open("figures/prog02/2_3_rbgs_k1_error.txt",
            ios::out | ios::trunc);
      k3file.open("figures/prog02/2_3_rbgs_k3_error.txt",
            ios::out | ios::trunc);
      k6file.open("figures/prog02/2_3_rbgs_k6_error.txt",
            ios::out | ios::trunc);

      // run methods and write to files
      if ( k1file.is_open() && k3file.is_open() && k6file.is_open() ) {

         rbgauss_seidel_write_errors(A,f,v1,100, k1file);
         rbgauss_seidel_write_errors(A,f,v3,100, k3file);
         rbgauss_seidel_write_errors(A,f,v6,100, k6file);

         k1file.close();
         k3file.close();
         k6file.close();
      }
      
      else {
         cerr << "file output error" << endl;
         exit(-1);
      }
   }

   else  {
      cerr << "method = " << method << " isn't implemented" << endl;
      exit(-1);
   }
}


int main(void) {
   
   make_2_3_plot_data(1);
   make_2_3_plot_data(2);
   make_2_3_plot_data(3);

   return 1;
}
