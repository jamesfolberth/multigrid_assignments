// prog03.cpp
//
// James Folberth
// Spring 2015

#include <cmath>
#include <iostream>
#include <fstream>
#include <functional>
#include <string>
#include <valarray>

#include "matrix_crs.hpp"
#include "model_problems.hpp"
#include "utils.hpp"
#include "multigrid.hpp"
using namespace std;

void table_4_1(int mode) {
   // {{{ 
   unsigned L, Lx, Ly; 
   string filename;
   switch (mode) {
      case 0:
         Lx = 4;
         Ly = 4;
         L = 2;
         filename = "figures/prog03/table_4_1_16.txt";
         break;

      case 1:
         Lx = 5;
         Ly = 5;
         L = 3;
         filename = "figures/prog03/table_4_1_32.txt";
         break;

      case 2:
         Lx = 6;
         Ly = 6;
         L = 4;
         filename = "figures/prog03/table_4_1_64.txt";
         break;

      case 3:
         Lx = 7;
         Ly = 7;
         L = 5;
         filename = "figures/prog03/table_4_1_128.txt";
         break;

      default:
         cerr << "error: prog03.cpp:table_4_1: bad mode" << endl;
         exit(-1);
   }
 
   unsigned n = (pow(2,Lx)-1)*(pow(2,Ly)-1);
   double sigma = 0., resid_nrm, err_nrm, resid_nrm_prev, err_nrm_prev;
   valarray<double> f(0.,n), u(0.,n), v0(0.,n), v(0.,n),
      resid(0.,n), err(0.,n);

   unsigned nx = pow(2,Lx)-1, ny = pow(2,Ly)-1;

   // Equation 4.8 of the text
   sigma = 0.;
   v0 = rand_vec<double>(n,-1.,1.);
   double xij, yij;
   for (unsigned i=0; i < ny; ++i) {
      for (unsigned j=0; j < nx; ++j) {
         xij = double(j+1.)/double(nx+1.);
         yij = double(i+1.)/double(ny+1.);
         u[ny*i+j] = (xij*xij - pow(xij,4.))*(pow(yij,4.)-yij*yij);
         f[ny*i+j] = 2.*((1.-6.*xij*xij)*yij*yij*(1.-yij*yij) 
               + (1.-6.*yij*yij)*xij*xij*(1.-xij*xij));
      }
   }

   // open file for writing
   fstream file;
   file.open(filename, ios::out | ios::trunc);
   if ( file.is_open() ) {
  
   // make a function to build A that depends only on the grid level
   // (so fix sigma in the beginning)
   // [=] pass by copy
   auto build_A = [=](unsigned _Lx, unsigned _Ly) -> matrix_crs<double>
      {return model_problem_2d<double>(_Lx,_Ly,sigma);};
   auto build_P = [](unsigned _Lx, unsigned _Ly) -> matrix_crs<double>
      {return operator_2d_interp_lin<double>(_Lx,_Ly);};
   auto build_R = [](unsigned _Lx, unsigned _Ly) -> matrix_crs<double>
      //{return operator_2d_restrict_inj<double>(_Lx,_Ly);};
      {return operator_2d_restrict_full<double>(_Lx,_Ly);};

   auto smoother = [](const matrix_crs<double>& _A, const valarray<double>& _f,
         valarray<double>& _v, unsigned _num_itr)
      //{wjacobi_ip<double>(_A,_f,_v,_num_itr);};
      //{gauss_seidel_ip<double>(_A,_f,_v,_num_itr);};
      {rbgauss_seidel_ip<double>(_A,_f,_v,_num_itr);};


   // build levels
   vector<level<double>> levels = build_levels_2d<double>(build_A, f, build_P,
         build_R, smoother, L, Lx, Ly, v0);

   // do V cycle
   file << "n = " << nx << " points in each direction" << endl;
   file << "Matrix A is " << levels[0].A.m << " x " << levels[0].A.n 
        << " with " << levels[0].A.val.size() << " nonzero elements" << endl;
   
   resid = levels[0].f - levels[0].A*levels[0].v;
   err = u-levels[0].v;
   resid_nrm_prev = dl2norm(resid, nx, ny);
   err_nrm_prev = dl2norm(err, nx, ny);
   
   file << "||r||_h = " << _PRINT_VECTOR_FORMAT_ << resid_nrm_prev
        << "  ||e||_h = " << _PRINT_VECTOR_FORMAT_ << err_nrm_prev << endl;
 
   for (unsigned i = 0; i < 15; ++i) {
      vcycle(levels, levels.begin(), 2, 1);
      //mucycle(levels, levels.begin(), 2, 1, 2);

      resid = levels[0].f - levels[0].A*levels[0].v;
      err = u-levels[0].v;
      resid_nrm = dl2norm(resid, nx, ny);
      err_nrm = dl2norm(err, nx, ny);

      file << "||r||_h = " << _PRINT_VECTOR_FORMAT_ << resid_nrm
           << "  ratio = " << _PRINT_VECTOR_FORMAT_ 
           << resid_nrm/resid_nrm_prev
           << "  ||e||_h = " << _PRINT_VECTOR_FORMAT_ << err_nrm
           << "  ratio = " << _PRINT_VECTOR_FORMAT_ 
           << err_nrm/err_nrm_prev << "  at i = " << i << endl;

      resid_nrm_prev = resid_nrm;
      err_nrm_prev = err_nrm;
   }

   // close file
   file.close();
   } // if ( file.is_open() )
   else {
      cerr << "error: prog03.cpp:table_4_1: file IO error" << endl;
      exit(-1);
   }

   // }}}
}

void table_4_2(void) {
   // {{{ 
   unsigned L = 7, Lx = 10, Ly = 10, n = (pow(2,Lx)-1)*(pow(2,Ly)-1);
   double sigma = 0., resid_nrm, resid_nrm_prev, resid_ratio_avg;
   valarray<double> f(0.,n), u(0.,n), v0(0.,n), v(0.,n),
      resid(0.,n), err(0.,n);

   unsigned nx = pow(2,Lx)-1, ny = pow(2,Ly)-1;
   string smoother_name, restrict_name;

   // Equation 4.8 of the text
   sigma = 0.;
   v0 = rand_vec<double>(n,-1.,1.);
   double xij, yij;
   for (unsigned i=0; i < ny; ++i) {
      for (unsigned j=0; j < nx; ++j) {
         xij = double(j+1.)/double(nx+1.);
         yij = double(i+1.)/double(ny+1.);
         u[ny*i+j] = (xij*xij - pow(xij,4.))*(pow(yij,4.)-yij*yij);
         f[ny*i+j] = 2.*((1.-6.*xij*xij)*yij*yij*(1.-yij*yij) 
               + (1.-6.*yij*yij)*xij*xij*(1.-xij*xij));
      }
   }

   // make a function to build A that depends only on the grid level
   // (so fix sigma in the beginning)
   // [=] pass by copy
   auto build_A = [=](unsigned _Lx, unsigned _Ly) -> matrix_crs<double>
      {return model_problem_2d<double>(_Lx,_Ly,sigma);};
   auto build_P = [](unsigned _Lx, unsigned _Ly) -> matrix_crs<double>
      {return operator_2d_interp_lin<double>(_Lx,_Ly);};
   
   //restrict_name = "Lin. Inj.";
   restrict_name = "Lin. Full.";
   auto build_R = [](unsigned _Lx, unsigned _Ly) -> matrix_crs<double>
      //{return operator_2d_restrict_inj<double>(_Lx,_Ly);};
      {return operator_2d_restrict_full<double>(_Lx,_Ly);};

   smoother_name = "Jacobi";
   //smoother_name = "GS";
   //smoother_name = "RBGS";
   auto smoother = [](const matrix_crs<double>& _A, const valarray<double>& _f,
         valarray<double>& _v, unsigned _num_itr)
      {wjacobi_ip<double>(_A,_f,_v,_num_itr);};
      //{gauss_seidel_ip<double>(_A,_f,_v,_num_itr);};
      //{rbgauss_seidel_ip<double>(_A,_f,_v,_num_itr);};


   // build levels
   vector<level<double>> levels = build_levels_2d<double>(build_A, f, build_P,
         build_R, smoother, L, Lx, Ly, v0);

   unsigned num_cycles = 5;
   unsigned nu1, nu2;
   
   resid = levels[0].f - levels[0].A*levels[0].v;
   resid_nrm_prev = dl2norm(resid, nx, ny);

   nu1 = 1; nu2 = 0;
   resid_ratio_avg = 1.;
   for (unsigned i = 0; i < num_cycles; ++i) {
      vcycle(levels, levels.begin(), nu1, nu2);

      resid = levels[0].f - levels[0].A*levels[0].v;
      resid_nrm = dl2norm(resid, nx, ny);

      resid_ratio_avg *= resid_nrm / resid_nrm_prev;
      resid_nrm_prev = resid_nrm;
      
   }

   cout << "V(" << nu1 << "," << nu2 << ") " << smoother_name << " "
        << restrict_name << " resid_ratio_avg = " << _PRINT_VECTOR_FORMAT_ 
        << pow(resid_ratio_avg, 1./static_cast<double>(num_cycles)) << endl;


   levels[0].v = v0;
   resid = levels[0].f - levels[0].A*levels[0].v;
   resid_nrm_prev = dl2norm(resid, nx, ny);


   nu1 = 1; nu2 = 1;
   resid_ratio_avg = 1.;
   for (unsigned i = 0; i < num_cycles; ++i) {
      vcycle(levels, levels.begin(), nu1, nu2);

      resid = levels[0].f - levels[0].A*levels[0].v;
      resid_nrm = dl2norm(resid, nx, ny);

      resid_ratio_avg *= resid_nrm / resid_nrm_prev;
      resid_nrm_prev = resid_nrm;
      
   }

   cout << "V(" << nu1 << "," << nu2 << ") " << smoother_name << " "
        << restrict_name << " resid_ratio_avg = " << _PRINT_VECTOR_FORMAT_ 
        << pow(resid_ratio_avg, 1./static_cast<double>(num_cycles)) << endl;


   levels[0].v = v0;
   resid = levels[0].f - levels[0].A*levels[0].v;
   resid_nrm_prev = dl2norm(resid, nx, ny);

   nu1 = 2; nu2 = 1;
   resid_ratio_avg = 1.;
   for (unsigned i = 0; i < num_cycles; ++i) {
      vcycle(levels, levels.begin(), nu1, nu2);

      resid = levels[0].f - levels[0].A*levels[0].v;
      resid_nrm = dl2norm(resid, nx, ny);

      resid_ratio_avg *= resid_nrm / resid_nrm_prev;
      resid_nrm_prev = resid_nrm;
      
   }

   cout << "V(" << nu1 << "," << nu2 << ") " << smoother_name << " "
        << restrict_name << " resid_ratio_avg = " << _PRINT_VECTOR_FORMAT_ 
        << pow(resid_ratio_avg, 1./static_cast<double>(num_cycles)) << endl;
 
   // }}}
}

void table_4_2_1d(void) {
   // {{{ 
   unsigned L = 7, Lx = 10, n = (pow(2,Lx)-1);
   double sigma = 0., resid_nrm, resid_nrm_prev, resid_ratio_avg;
   valarray<double> f(0.,n), u(0.,n), v0(0.,n), v(0.,n),
      resid(0.,n), err(0.,n);

   unsigned nx = pow(2,Lx)-1;
   string smoother_name, restrict_name;
   
   double xi = 0., C = 2., k = 3.;

   int mode = 1;
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
         sigma = 1.;
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
   // [=] pass by copy
   auto build_A = [=](unsigned _Lx) -> matrix_crs<double>
      {return model_problem_1d<double>(_Lx,sigma);};
   auto build_P = [](unsigned _Lx) -> matrix_crs<double>
      {return operator_1d_interp_lin<double>(_Lx);};
   
   //restrict_name = "Lin. Inj.";
   restrict_name = "Lin. Full.";
   auto build_R = [](unsigned _Lx) -> matrix_crs<double>
      //{return operator_1d_restrict_inj<double>(_Lx);};
      {return operator_1d_restrict_full<double>(_Lx);};

   //smoother_name = "Jacobi";
   //smoother_name = "GS";
   smoother_name = "RBGS";
   auto smoother = [](const matrix_crs<double>& _A, const valarray<double>& _f,
         valarray<double>& _v, unsigned _num_itr)
      //{wjacobi_ip<double>(_A,_f,_v,_num_itr);};
      //{gauss_seidel_ip<double>(_A,_f,_v,_num_itr);};
      {rbgauss_seidel_ip<double>(_A,_f,_v,_num_itr);};


   // build levels
   vector<level<double>> levels = build_levels_1d<double>(build_A, f, build_P,
         build_R, smoother, L, Lx, v0, 0);

   unsigned num_cycles = 10;
   unsigned nu1, nu2;
   
   resid = levels[0].f - levels[0].A*levels[0].v;
   resid_nrm_prev = dl2norm(resid, nx);

   nu1 = 1; nu2 = 0;
   resid_ratio_avg = 1.;
   for (unsigned i = 0; i < num_cycles; ++i) {
      vcycle(levels, levels.begin(), nu1, nu2);

      resid = levels[0].f - levels[0].A*levels[0].v;
      resid_nrm = dl2norm(resid, nx);

      resid_ratio_avg *= resid_nrm / resid_nrm_prev;
      resid_nrm_prev = resid_nrm;
      
   }

   cout << "V(" << nu1 << "," << nu2 << ") " << smoother_name << " "
        << restrict_name << " resid_ratio_avg = " << _PRINT_VECTOR_FORMAT_ 
        << pow(resid_ratio_avg, 1./static_cast<double>(num_cycles)) << endl;


   levels[0].v = v0;
   resid = levels[0].f - levels[0].A*levels[0].v;
   resid_nrm_prev = dl2norm(resid, nx);


   nu1 = 1; nu2 = 1;
   resid_ratio_avg = 1.;
   for (unsigned i = 0; i < num_cycles; ++i) {
      vcycle(levels, levels.begin(), nu1, nu2);

      resid = levels[0].f - levels[0].A*levels[0].v;
      resid_nrm = dl2norm(resid, nx);

      resid_ratio_avg *= resid_nrm / resid_nrm_prev;
      resid_nrm_prev = resid_nrm;
      
   }

   cout << "V(" << nu1 << "," << nu2 << ") " << smoother_name << " "
        << restrict_name << " resid_ratio_avg = " << _PRINT_VECTOR_FORMAT_ 
        << pow(resid_ratio_avg, 1./static_cast<double>(num_cycles)) << endl;


   levels[0].v = v0;
   resid = levels[0].f - levels[0].A*levels[0].v;
   resid_nrm_prev = dl2norm(resid, nx);

   nu1 = 2; nu2 = 1;
   resid_ratio_avg = 1.;
   for (unsigned i = 0; i < num_cycles; ++i) {
      vcycle(levels, levels.begin(), nu1, nu2);

      resid = levels[0].f - levels[0].A*levels[0].v;
      resid_nrm = dl2norm(resid, nx);

      resid_ratio_avg *= resid_nrm / resid_nrm_prev;
      resid_nrm_prev = resid_nrm;
      
   }

   cout << "V(" << nu1 << "," << nu2 << ") " << smoother_name << " "
        << restrict_name << " resid_ratio_avg = " << _PRINT_VECTOR_FORMAT_ 
        << pow(resid_ratio_avg, 1./static_cast<double>(num_cycles)) << endl;
 
   // }}}
}

void mg_1d_vcycle(void) {
   // {{{ 
   unsigned L = 7, Lx = 10, n = pow(2,Lx)-1;
   double sigma = 0., resid_nrm, err_nrm, resid_nrm_prev, resid_ratio_avg;
   valarray<double> f(0.,n), u(0.,n), v0(0.,n), v(0.,n),
      resid(0.,n), err(0.,n);

   double xi = 0.,C = 2., k = 3.;

   int mode = 1;
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
   vector<level<double>> levels = build_levels_1d<double>(build_A, f, build_P,
         build_R, smoother, L, Lx, v0, 0);


   // Print to a file
   fstream file;
   file.open("figures/prog03/mg_1d_vcycle_J.txt", ios::out | ios::trunc);
   //file.open("figures/prog03/mg_1d_vcycle_GS.txt", ios::out | ios::trunc);
   //file.open("figures/prog03/mg_1d_vcycle_RBGS.txt", ios::out | ios::trunc);
   if ( file.is_open() ) {

   // do V cycle
   file << "Matrix A is " << levels[0].A.m << " x " << levels[0].A.n << endl;
   
   resid = levels[0].f - levels[0].A*levels[0].v;
   err = u-levels[0].v;
   resid_nrm_prev = dl2norm(resid,n);
   err_nrm = dl2norm(err,n);
   
   file << "||r||_h = " << _PRINT_VECTOR_FORMAT_ << resid_nrm_prev
        << "  ||e||_h = " << _PRINT_VECTOR_FORMAT_ << err_nrm << endl;

   resid_ratio_avg = 1.;
   unsigned num_cycles = 10;
   for (unsigned i = 0; i < num_cycles; ++i) {
      vcycle(levels, levels.begin(), 2, 1);

      resid = levels[0].f - levels[0].A*levels[0].v;
      err = u-levels[0].v;
      resid_nrm = dl2norm(resid,n);
      err_nrm = dl2norm(err,n);

      file << "||r||_h = " << _PRINT_VECTOR_FORMAT_ << resid_nrm
           << "  ratio = " << _PRINT_VECTOR_FORMAT_ 
           << resid_nrm/resid_nrm_prev
           << "  ||e||_h = " << _PRINT_VECTOR_FORMAT_ << err_nrm
           << "  at i = " << i << endl;

      resid_ratio_avg *= resid_nrm / resid_nrm_prev;
      resid_nrm_prev = resid_nrm;
   }

   file << "Average reduction fator = " 
        << pow(resid_ratio_avg, 1./static_cast<double>(num_cycles)) << endl;

   }
   else { 
      cerr << "error: prog03.cpp:mg_1d_vcycle: file IO error" << endl;
      exit(-1);
   }

   // }}}
}

int main(void) {
   
   table_4_1(0);
   table_4_1(1);
   table_4_1(2);
   table_4_1(3);

   //table_4_2();// need to change function and run 6 times
   //table_4_2_1d();// need to change function and run 6 times
   
   mg_1d_vcycle();// need to change function and run 3 times

   return 0;
}
