#include <cmath>
#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <random>

#include <ext/numeric>
using __gnu_cxx::power;
#include <vector>
#include <utility>

#include <boost/python.hpp>
#include <boost/math/constants/constants.hpp>

#include "numpy_bind.hh"


using namespace std;
using namespace boost;

std::mt19937 gen{42};
const double pi = boost::math::constants::pi<double>();

void get_alpha(python::object oralpha, python::object oialpha,
               python::object of_i, python::object of_range, double fmin)
{
    multi_array_ref<double, 2> ralpha = get_array<double,2>(oralpha);
    multi_array_ref<double, 2> ialpha = get_array<double,2>(oialpha);

    multi_array_ref<double, 1> f_i = get_array<double,1>(of_i);
    multi_array_ref<double, 1> f_range = get_array<double,1>(of_range);

    int i;
    for (i = 0; i < int(ralpha.shape()[0]); ++i)
    {
        for (size_t j = 0; j < ralpha.shape()[1]; ++j)
        {
	 
	  ralpha[i][j] = ralpha[int(f_i[int(-f_range[i] - fmin)])]
	    [int(f_i[int(-f_range[j] - fmin)])];
	  ialpha[i][j] = -ialpha[int(f_i[int(-f_range[i] - fmin)])]
	    [int(f_i[int(-f_range[j] - fmin)])];
        }
    }
}

void get_v_pos(int i, int j, multi_array_ref<double, 2>& eta,
               double& vx, double& vy, double delta)
{
    int jp = (j < int(eta.shape()[1] - 1)) ? (j + 1) : 0;
    int jm = (j > 0) ? j - 1 : eta.shape()[1] - 1;
  
    vx = (eta[i][jp] - eta[i][jm])/(delta*2.);


    int ip = (i < int(eta.shape()[0] - 1)) ? (i + 1) : 0;
    int im = (i > 0) ? i - 1 : eta.shape()[0] - 1;
    vy = -(eta[ip][j] - eta[im][j])/(delta*2.);
}

void get_vel(python::object oeta, python::object ovx, python::object ovy, double delta)
{
    multi_array_ref<double, 2> eta = get_array<double,2>(oeta);
    multi_array_ref<double, 2> vx = get_array<double,2>(ovx);
    multi_array_ref<double, 2> vy = get_array<double,2>(ovy);

    int i;
    for (i = 0; i < int(eta.shape()[0]); ++i)
    {
        for (int j = 0; j < int(eta.shape()[1]); ++j)
        {
            get_v_pos(i, j, eta, vx[i][j], vy[i][j], delta);
        }
    }
}




void vorticity(python::object os, python::object oscos, python::object ossin,
               python::object ovx, python::object ovy, double delta)
{
    multi_array_ref<double, 2> s = get_array<double,2>(os);
    multi_array_ref<double, 2> scos = get_array<double,2>(oscos);
    multi_array_ref<double, 2> ssin = get_array<double,2>(ossin);
    multi_array_ref<double, 2> vx = get_array<double,2>(ovx);
    multi_array_ref<double, 2> vy = get_array<double,2>(ovy);

    int i;
    
    for (i = 0; i < int(s.shape()[0]); ++i)
    {
        for (int j = 0; j < int(s.shape()[1]); ++j)
        {
            int ip = (i < int(s.shape()[0] - 1)) ? (i + 1) : 0;
            int im = (i > 0) ? i - 1 : s.shape()[0] - 1;

            int jp = (j < int(s.shape()[1] - 1)) ? (j + 1) : 0;
            int jm = (j > 0) ? j - 1 : s.shape()[1] - 1;

            double sxy, syx, sxx, syy;
            sxx = (vx[ip][j]-vx[im][j])/(2.*delta);
            syy = (vy[i][jp]-vy[i][jm])/(2.*delta);
            sxy = (vx[i][jp]-vx[i][jm])/(2.*delta);
            syx = (vy[ip][j]-vy[im][j])/(2.*delta);

            
            s[i][j] = -(sxy - syx)/2.;
            scos[i][j] = (syx + sxy)/2.;
            ssin[i][j] = (sxx - syy)/2.;
            
            
              
        }
    }   
}

void move_bacteria(int N_bac, python::object oux, python::object ouy,
                   python::object oomega,  python::object oscos,
                   python::object ossin, python::object ox,
                   python::object odx,  python::object ody, python::object odtheta,
                   double delta, int N, double alpha, double us,
                   double ug)
{

    multi_array_ref<double, 1> x = get_array<double,1>(ox);
    multi_array_ref<double, 1> dx = get_array<double,1>(odx);
    multi_array_ref<double, 1> dy = get_array<double,1>(ody);
    multi_array_ref<double, 1> dtheta = get_array<double,1>(odtheta);
    multi_array_ref<double, 2> ux = get_array<double,2>(oux);
    multi_array_ref<double, 2> uy = get_array<double,2>(ouy);
    multi_array_ref<double, 2> omega = get_array<double,2>(oomega);
    multi_array_ref<double, 2> ssin = get_array<double,2>(ossin);
    multi_array_ref<double, 2> scos = get_array<double,2>(oscos);
    
    
     for (int i = 0; i < N_bac; ++i)
     {
         double x1 = x[i];
         double x2 = x[N_bac + i];
         double theta = x[2*N_bac + i];

         int xn = fmod(int(x1/delta), N);
         int yn = fmod(int(x2/delta), N);

         double ux_pos = ux[xn][yn];
         double uy_pos = uy[xn][yn];
         
         dtheta[i] = omega[xn][yn] + alpha*(scos[xn][yn]*cos(2*theta)  - ssin[xn][yn]*sin(2*theta));
         dx[i] = ux_pos + us*cos(theta);
         dy[i] = uy_pos + us*sin(theta) + ug;
         
     }
     

}

BOOST_PYTHON_MODULE(libfluid)
{
    using namespace boost::python;
    def("get_alpha", get_alpha);
    def("get_vel", get_vel);
    def("vorticity", vorticity);
    def("move_bacteria", move_bacteria);
    
};
