// prog04.cpp
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

// 1D
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
   vector<level<double>> levels = build_levels_1d<double>(build_A, f, build_P,
         build_R, smoother, L, Lx, v0, 0);

   // do V cycle
   cout << "Matrix A is " << levels[0].A.m << " x " << levels[0].A.n << endl;
   
   resid = levels[0].f - levels[0].A*levels[0].v;
   err = u-levels[0].v;
   resid_nrm_prev = dl2norm(resid,n);
   err_nrm = dl2norm(err,n);
   
   cout << "||r||_h = " << _PRINT_VECTOR_FORMAT_ << resid_nrm_prev
        << "  ||e||_h = " << _PRINT_VECTOR_FORMAT_ << err_nrm << endl;

   resid_ratio_avg = 1.;
   unsigned num_cycles = 10;
   for (unsigned i = 0; i < num_cycles; ++i) {
      vcycle(levels, levels.begin(), 2, 1);

      resid = levels[0].f - levels[0].A*levels[0].v;
      err = u-levels[0].v;
      resid_nrm = dl2norm(resid,n);
      err_nrm = dl2norm(err,n);

      cout << "||r||_h = " << _PRINT_VECTOR_FORMAT_ << resid_nrm
           << "  ratio = " << _PRINT_VECTOR_FORMAT_ 
           << resid_nrm/resid_nrm_prev
           << "  ||e||_h = " << _PRINT_VECTOR_FORMAT_ << err_nrm
           << "  at i = " << i << endl;

      resid_ratio_avg *= resid_nrm / resid_nrm_prev;
      resid_nrm_prev = resid_nrm;
   }

   cout << "Average reduction fator = " 
        << pow(resid_ratio_avg, 1./static_cast<double>(num_cycles)) << endl;

   // }}}
}

void mg_1d_fmg(void) {
   // {{{ 
   unsigned L = 7, Lx = 10, n = pow(2,Lx)-1;
   double sigma = 0., resid_nrm, err_nrm;
   valarray<double> f(0.,n), u(0.,n), u_coarse(0.,pow(2,Lx-L+1)-1),
      v0(0.,n), v(0.,n), resid(0.,n), err(0.,n);

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
   // Note: this will create a random guess on the finest grid, but we won't 
   // use it, since we're doing FMG.
   vector<level<double>> levels = build_levels_1d<double>(build_A, build_P,
         build_R, smoother, L, Lx);


   // Assign initial condition on coarsest grid and RHS vectors on all grids
   double xi = 0.,C = 2., k = 3.;
   unsigned n_temp;
   int mode = 1;
   switch (mode) {
      case 0:
         // zero RHS
   
         // true solution
         u *= 0.;
         u_coarse *= 0.;
         
         // Initial guess on coarse grid
         cout << levels[L-1].v.size() << endl;
         cout << pow(2,Lx-L+1)-1 << endl;
         levels[L-1].v = rand_vec<double>(pow(2,Lx-L+1)-1,-1.,1.);

         // Set RHS for all grids
         for (unsigned l=0; l < L; ++l) {
            n_temp = pow(2,Lx-l)-1;

            for (unsigned i=0; i < n_temp; ++i) {
               levels[l].f[i] = 0.;
            }
         }
         break;

      case 1:
         // sin(pi*k*x) RHS

         // True solution
         for (unsigned i=0; i < n; ++i) {
            u[i] = C/(pow(_PI_*k,2.)+sigma)*sin(k*_PI_*double(i+1-0)/double(n+1));
         }

         n_temp = pow(2,Lx-L+1)-1;
         for (unsigned i=0; i < n_temp; ++i) {
            u_coarse[i] = C/(pow(_PI_*k,2.)+sigma)*sin(k*_PI_*double(i+1-0)
                  /double(n_temp+1));
         }

         // Initial guess on coarse grid
         levels[L-1].v = rand_vec<double>(pow(2,Lx-L+1)-1,-1.,1.);

         // Set RHS for all grids
         for (unsigned l=0; l < L; ++l) {
            n_temp = pow(2,Lx-l)-1;

            for (unsigned i=0; i < n_temp; ++i) {
               levels[l].f[i] = C*sin(k*_PI_*double(i+1-0)/double(n_temp+1));
            }

            // neet to account for Dirichlet BCs
            levels[l].f[0] += pow(n_temp+1,2)*0.;
            levels[l].f[n_temp-1] += pow(n_temp+1,2)*0.;
         }

         break;

      case 2:
         // sin(pi*x)*exp(-x) RHS
         
         // True solution
         for (unsigned i=0; i < n; ++i) {
            xi = double(i+1-0)/double(n+1);
            u[i] = sin(_PI_*xi)*exp(-xi);
         }

         n_temp = pow(2,Lx-L+1)-1;
         for (unsigned i=0; i < n_temp; ++i) {
            xi = double(i+1-0)/double(n_temp+1);
            u_coarse[i] = sin(_PI_*xi)*exp(-xi);
         }

         // Initial guess on coarse grid
         levels[L-1].v = rand_vec<double>(n_temp,-1.,1.);

         // Set RHS for all grids
         for (unsigned l=0; l < L; ++l) {
            n_temp = pow(2,Lx-l)-1;

            for (unsigned i=0; i < n_temp; ++i) {
               xi = double(i+1-0)/double(n_temp+1);
               levels[l].f[i] = exp(-xi)*((-1+_PI_*_PI_+sigma)*sin(_PI_*xi) 
                     + 2.*_PI_*cos(_PI_*xi));
 
            }

            // neet to account for Dirichlet BCs
            levels[l].f[0] += pow(n_temp+1,2)*0.;
            levels[l].f[n_temp-1] += pow(n_temp+1,2)*0.;
         }

         break;
 
      default:
         cerr << "error: test.cpp:test_mg_1d_fmg: bad mode" << endl;
         exit(-1);
   }

  
   // do FMG
   cout << "Top-level A is " << levels[0].A.m << " x " << levels[0].A.n 
        << " with " << levels[0].A.val.size() << " nonzero elements" << endl;
 
   resid = levels[L-1].f - levels[L-1].A*levels[L-1].v;
   err = u_coarse-levels[L-1].v;
   resid_nrm = dl2norm(resid,pow(2,Lx-L+1)-1);
   err_nrm = dl2norm(err,pow(2,Lx-L+1)-1);
   
   cout << "Pre FMG:  ||r||_h = " << _PRINT_VECTOR_FORMAT_ << resid_nrm
        << "  ||e||_h = " << _PRINT_VECTOR_FORMAT_ << err_nrm << endl;

   fmg(levels, 2, 1, 1);

   resid = levels[0].f - levels[0].A*levels[0].v;
   err = u-levels[0].v;
   resid_nrm = dl2norm(resid,n);
   err_nrm = dl2norm(err,n);
   
   cout << "Post FMG: ||r||_h = " << _PRINT_VECTOR_FORMAT_ << resid_nrm
        << "  ||e||_h = " << _PRINT_VECTOR_FORMAT_ << err_nrm << endl;


   // }}}
}

// 2D
void mg_2d_vcycle(void) {
   // {{{ 
   unsigned L = 7, Lx = 10, Ly = 10, n = (pow(2,Lx)-1)*(pow(2,Ly)-1);
   double sigma = 2., resid_nrm, err_nrm, resid_nrm_prev;
   valarray<double> f(0.,n), u(0.,n), v0(0.,n), v(0.,n),
      resid(0.,n), err(0.,n);

   unsigned nx = pow(2,Lx)-1, ny = pow(2,Ly)-1;

   int mode = 1;
   switch (mode) {
      case 0:
         // zero RHS
         sigma = 0.;
         v0 = rand_vec<double>(n,-1.,1.);
         f *= 0.; u *= 0.;
         break;

      case 1:
      {
         // sin(pi*kx*x)*sin(pi*ky*y) RHS
         unsigned kx = 3, ky = 10;
         double C = 3.;

         v0 = rand_vec<double>(n,-1.,1.);
         double xij, yij;
         for (unsigned i=0; i < ny; ++i) {
            for (unsigned j=0; j < nx; ++j) {
               xij = double(j+1.)/double(nx+1.);
               yij = double(i+1.)/double(ny+1.);
               u[ny*i+j] = C/(pow(_PI_*kx,2.)+pow(_PI_*ky,2.)+sigma)
                  *sin(_PI_*kx*xij)*sin(_PI_*ky*yij);
               f[ny*i+j] = C*sin(_PI_*kx*xij)*sin(_PI_*ky*yij);
            }
         }
         break;
      }

      case 2:
      {
         // Equation 4.8 from the text
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
         break;
      }
      default:
         cerr << "error: test.cpp:test_mg_2d_vcycle: bad mode" << endl;
         exit(-1);
   }

  
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
   cout << "Matrix A is " << levels[0].A.m << " x " << levels[0].A.n 
        << " with " << levels[0].A.val.size() << " nonzero elements" << endl;
   
   resid = levels[0].f - levels[0].A*levels[0].v;
   err = u-levels[0].v;
   resid_nrm_prev = dl2norm(resid, nx, ny);
   err_nrm = dl2norm(err, nx, ny);
   
   cout << "||r||_h = " << _PRINT_VECTOR_FORMAT_ << resid_nrm_prev
        << "  ||e||_h = " << _PRINT_VECTOR_FORMAT_ << err_nrm << endl;
 
   for (unsigned i = 0; i < 10; ++i) {
      vcycle(levels, levels.begin(), 2, 1);
      //mucycle(levels, levels.begin(), 2, 1, 2);

      resid = levels[0].f - levels[0].A*levels[0].v;
      err = u-levels[0].v;
      resid_nrm = dl2norm(resid, nx, ny);
      //resid_nrm = norm(resid,0);
      err_nrm = dl2norm(err, nx, ny);
      //err_nrm = norm(err,0);

      cout << "||r||_h = " << _PRINT_VECTOR_FORMAT_ << resid_nrm
           << "  ||e||_h = " << _PRINT_VECTOR_FORMAT_ << err_nrm
           << "  ratio = " << _PRINT_VECTOR_FORMAT_ 
           << resid_nrm/resid_nrm_prev << "  at i = " << i << endl;

      resid_nrm_prev = resid_nrm;
   }

   // }}}
}

void mg_2d_fmg(void) {
   // {{{ 
   unsigned L = 7, Lx = 10, Ly = 10, n = (pow(2,Lx)-1)*(pow(2,Ly)-1);
   double sigma = 2., resid_nrm, err_nrm;
   valarray<double> f(0.,n), u(0.,n),
      u_coarse(0.,(pow(2,Lx-L+1)-1)*(pow(2,Ly-L+1)-1)),
      v0(0.,n), v(0.,n), resid(0.,n), err(0.,n);

   unsigned nx = pow(2,Lx)-1, ny = pow(2,Ly)-1, nx_temp, ny_temp;


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
      //{wjacobi_ip<double>(_A,_f,_v,_num_itr, 4./5.);};
      //{gauss_seidel_ip<double>(_A,_f,_v,_num_itr);};
      {rbgauss_seidel_ip<double>(_A,_f,_v,_num_itr);};


   // build levels
   vector<level<double>> levels = build_levels_2d<double>(build_A, build_P,
         build_R, smoother, L, Lx, Ly);

   int mode = 1;
   switch (mode) {
      case 0:
         // zero RHS

         // True solution
         u *= 0.;
         u_coarse *= 0.;

         // Initial guess on coarse grid
         levels[L-1].v = rand_vec<double>((pow(2,Lx-L+1)-1)*(pow(2,Ly-L+1)-1),-1.,1.);

         // Set RHS for all grids
         for (unsigned l=0; l < L; ++l) {
            nx_temp = pow(2,Lx-l)-1;
            ny_temp = pow(2,Ly-l)-1;

            for (unsigned i=0; i < ny_temp; ++i) {
               for (unsigned j=0; j < nx_temp; ++j) {
                  levels[l].f[ny_temp*i+j] = 0.;
               }
            }
         }


         break;

      case 1:
      {
         // sin(pi*k*x)*sin(pi*l*y) RHS
         unsigned kx = 3, ky = 10;
         double C = 3.;
         double xij, yij;

         // True solution
         for (unsigned i=0; i < ny; ++i) {
            for (unsigned j=0; j < nx; ++j) {
               xij = double(j+1.)/double(nx+1.);
               yij = double(i+1.)/double(ny+1.);
               u[ny*i+j] = C/(pow(_PI_*kx,2.)+pow(_PI_*ky,2.)+sigma)
                  *sin(_PI_*kx*xij)*sin(_PI_*ky*yij);
 
            }
         }

         nx_temp = pow(2,Lx-L+1)-1;
         ny_temp = pow(2,Ly-L+1)-1;
         for (unsigned i=0; i < ny_temp; ++i) {
            for (unsigned j=0; j < nx_temp; ++j) {
               xij = double(j+1.)/double(nx_temp+1.);
               yij = double(i+1.)/double(ny_temp+1.);
               u_coarse[ny*i+j] = C/(pow(_PI_*kx,2.)+pow(_PI_*ky,2.)+sigma)
                  *sin(_PI_*kx*xij)*sin(_PI_*ky*yij);
 
            }
         }

         // Initial guess on coarse grid
         levels[L-1].v = rand_vec<double>(nx_temp*ny_temp,-1.,1.);

         // Set RHS for all grids
         for (unsigned l=0; l < L; ++l) {
            nx_temp = pow(2,Lx-l)-1;
            ny_temp = pow(2,Ly-l)-1;

            for (unsigned i=0; i < ny_temp; ++i) {
               for (unsigned j=0; j < nx_temp; ++j) {
                  xij = double(j+1.)/double(nx_temp+1.);
                  yij = double(i+1.)/double(ny_temp+1.);
    
                  levels[l].f[ny_temp*i+j] = 
                     C*sin(_PI_*kx*xij)*sin(_PI_*ky*yij);
               }
            }
         }

         break;
      }
 

      case 2:
      {
         // Equation 4.8 from the text
         double xij, yij;

         // True solution
         for (unsigned i=0; i < ny; ++i) {
            for (unsigned j=0; j < nx; ++j) {
               xij = double(j+1.)/double(nx+1.);
               yij = double(i+1.)/double(ny+1.);
               u[ny*i+j] = (xij*xij - pow(xij,4.))*(pow(yij,4.)-yij*yij);
 
            }
         }

         nx_temp = pow(2,Lx-L+1)-1;
         ny_temp = pow(2,Ly-L+1)-1;
         for (unsigned i=0; i < ny_temp; ++i) {
            for (unsigned j=0; j < nx_temp; ++j) {
               xij = double(j+1.)/double(nx_temp+1.);
               yij = double(i+1.)/double(ny_temp+1.);
               u_coarse[ny*i+j] = (xij*xij - pow(xij,4.))*(pow(yij,4.)-yij*yij);
            }
         }

         // Initial guess on coarse grid
         levels[L-1].v = rand_vec<double>(nx_temp*ny_temp,-1.,1.);

         // Set RHS for all grids
         for (unsigned l=0; l < L; ++l) {
            nx_temp = pow(2,Lx-l)-1;
            ny_temp = pow(2,Ly-l)-1;

            for (unsigned i=0; i < ny_temp; ++i) {
               for (unsigned j=0; j < nx_temp; ++j) {
                  xij = double(j+1.)/double(nx_temp+1.);
                  yij = double(i+1.)/double(ny_temp+1.);
    
                  levels[l].f[ny_temp*i+j] = 
                     2.*((1.-6.*xij*xij)*yij*yij*(1.-yij*yij) 
                     + (1.-6.*yij*yij)*xij*xij*(1.-xij*xij));
               }
            }
         }

         break;
      }
      default:
         cerr << "error: test.cpp:test_mg_2d_fmg: bad mode" << endl;
         exit(-1);
   }

  
   // do FMG
   cout << "Top-level A is " << levels[0].A.m << " x " << levels[0].A.n 
        << " with " << levels[0].A.val.size() << " nonzero elements" << endl;
 
  
   resid = levels[L-1].f - levels[L-1].A*levels[L-1].v;
   err = u_coarse-levels[L-1].v;
   resid_nrm = dl2norm(resid, pow(2,Lx-L+1)-1, pow(2,Ly-L+1)-1);
   err_nrm = dl2norm(err, pow(2,Lx-L+1)-1, pow(2,Ly-L+1)-1);
   
   cout << "Pre FMG (coarse): ||r||_h = " << _PRINT_VECTOR_FORMAT_ << resid_nrm
        << "  ||e||_h = " << _PRINT_VECTOR_FORMAT_ << err_nrm << endl;

   fmg(levels, 2, 1, 1);

   resid = levels[0].f - levels[0].A*levels[0].v;
   err = u-levels[0].v;
   resid_nrm = dl2norm(resid, nx, ny);
   err_nrm = dl2norm(err, nx, ny);
   
   cout << "Post FMG (fine):  ||r||_h = " << _PRINT_VECTOR_FORMAT_ << resid_nrm
        << "  ||e||_h = " << _PRINT_VECTOR_FORMAT_ << err_nrm << endl;


   // }}}
}

int main(void) {

   //mg_1d_vcycle();
   //mg_1d_fmg();

   mg_2d_vcycle();
   mg_2d_fmg();

   return 0;
}
