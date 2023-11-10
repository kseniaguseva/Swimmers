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
    //#pragma omp parallel for default(shared) private(i) schedule(dynamic)
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
    //#pragma omp parallel for default(shared) private(i) schedule(dynamic)
    // cout << " c code" << "\n";
    // cout <<"shape " << int(eta.shape()[0]) << "  " << int(eta.shape()[1]) << "\n";
    for (i = 0; i < int(eta.shape()[0]); ++i)
    {
        for (int j = 0; j < int(eta.shape()[1]); ++j)
        {
            //cout << " "<< eta[i][j];
            get_v_pos(i, j, eta, vx[i][j], vy[i][j], delta);
            
        }
        //cout << "\n";
        
    }
    //cout << "c value " << eta[0][1] << "\n";
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
    
    //#pragma omp parallel for default(shared) private(i) schedule(dynamic)
    //cout << " c code \n";
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
            //cout << "  " << vx[i][j] <<" ";
            
              
        }
        //cout << "\n";
    }   

}

    
void get_interpolation(int Nx, int Ny, python::object oxf, python::object oyf,
                       double init_x, double init_y, python::object ox_range,
                       python::object oy_range, double dx, double dy,
                       python::object oMat, python::object oMat_new)
{
    multi_array_ref<double, 2> Mat = get_array<double,2>(oMat);
    multi_array_ref<double, 2> Mat_new = get_array<double,2>(oMat_new);
    multi_array_ref<double, 2> xf = get_array<double,2>(oxf);
    multi_array_ref<double, 2> yf = get_array<double,2>(oyf);
    multi_array_ref<double, 1> x_range = get_array<double,1>(ox_range);
    multi_array_ref<double, 1> y_range = get_array<double,1>(oy_range);

    int i, j;    

    for (i = 0; i < Nx; ++i)
    {
        for (j = 0; j < Ny; ++j)
        {
            double L = x_range[Nx-1] - x_range[0];
            double xp = xf[i][j];
            double yp = yf[i][j];
            
            xp = (xp <= L) ? xp : xp - L;
            xp = (xp >= 0) ? xp : L + xp;
            

            yp = (yp <= L) ? yp : yp - L;
            yp = (yp >= 0) ? yp : L + yp;

            
            int i1 = int(floor((xp - init_x)/dx));
            int j1 = int(floor((yp - init_y)/dy));
            
            i1 = (i1 >= 0) ? i1 : Nx-1;
            j1 = (j1 >= 0) ? j1 : Ny-1;
            
            int i2 = (i1 < int(Mat.shape()[0] - 1)) ? i1 + 1 : 0;
            int j2 = (j1 < int(Mat.shape()[1] - 1)) ? j1 + 1 : 0;            
            
            double x1 = x_range[i1];
            double x2 = x_range[i2];
            double y1 = y_range[j1];
            double y2 = y_range[j2];
            
            double c11 =  Mat[i1][j1];
            double c12 =  Mat[i1][j2];
            double c21 =  Mat[i2][j1];
            double c22 =  Mat[i2][j2];


            double delta = x_range[1] - x_range[0];
            //cout << "L " << Mat.shape()[0] <<"\n";
            // cout << "delta " << delta <<"\n";
            // cout << "L/delta " << L/delta <<"\n";
            // cout << "Nx " << Nx <<"\n";
            // cout << "dy " << dy <<"\n";

            if(x2 < xp){
                cout << "x2 " << x2 << "xp" << xp <<"\n";
            }
            
            
            // double delta_x1 = (x1 <= xp) ? xp - x1 : (L + xp) - x1;
            // double delta_y1 = (y1 <= yp) ? yp - y1 : (L + yp) - y1;

            // double delta_x2 = (x2 >= xp) ? x2 - xp : L + delta - xp;
            // double delta_y2 = (y2 >= yp) ? y2 - yp : L + delta - yp;

            double delta_x1 = xp - x1;
            double delta_y1 = yp - y1;

            double delta_x2 = x2 - xp;
            double delta_y2 = y2 - yp;

            
            //cout << "delta x2 " << delta_x2 <<"\n";
          
            
            Mat_new[i][j] = (c11*delta_x2*delta_y2 + c12*delta_x2*delta_y1
                             + c21*delta_x1*delta_y2 + c22*delta_x1*delta_y1)/(dx*dy);
          
            
        }
    }    
}


void get_concentration(double N_bac, python::object oc_last, python::object oC,
                         python::object ox,python::object oy,
                         double delta, int N)
{
    int i;

    multi_array_ref<double, 1> x = get_array<double, 1>(ox);
    multi_array_ref<double, 1> y = get_array<double, 1>(oy);
    
    multi_array_ref<double, 2> c_last = get_array<double, 2>(oc_last);
    multi_array_ref<double, 2> C = get_array<double, 2>(oC);
    for (i = 1; i < N_bac; ++i)
    {
        int xn, yn;
        xn = fmod((x[i]/delta),(N));
        yn = fmod((y[i]/delta),(N));
        c_last[0][i] = c_last[1][i];
        c_last[1][i] = C[xn][yn];
        
    }
    
}


void chemotaxis_all(python::object oc_last, python::object otheta, double Kd,
                      double alpha, double N_bac,
                      double sigma_rot, double tau_run, double dt, double phi0)
{
    int i;
    double Pr, event, ch, dp;
    double a, b;
    
    
    std::normal_distribution<> normal{0,1};
    std::uniform_real_distribution<> u{0,1};
    
    multi_array_ref<double, 2> c_last = get_array<double, 2>(oc_last);
    multi_array_ref<double, 1> theta = get_array<double, 1>(otheta);
    for (i = 1; i < N_bac; ++i)
    {        
            
        Pr = 0;
        
        // if(c_last[1] > c_last[0])
        // {
        dp = (c_last[i][1] - c_last[i][0])*(Kd/(pow((Kd + c_last[i][0]), 2)));
        ch = tau_run*exp(alpha*dp);
        // }
        // else
        // {
        //     ch = tau_run;
        // }

        
        Pr = (ch >1e-5)? dt/ch: 0;
        event = u(gen);
        if (event < Pr)
        {
            
            a = u(gen);
            //cout << "a" << a << "\n";
            a = (a > 0.5)? 1: -1;
            b = normal(gen)*sigma_rot;
            //cout << b << "\n";
            theta[i] = fmod(theta[i] + phi0*a + b, 2*pi);
            
        }   
    }
}



double get_chemotaxis(python::object ot_range, python::object oc_path,
                      double Kd, double T, double t, double tf)

{
    multi_array_ref<double, 1> t_range = get_array<double,1>(ot_range);
    multi_array_ref<double, 1> c_path = get_array<double,1>(oc_path);


    double dp;
    double dpdt_mean;
    dp = 0;
    dpdt_mean = 0;
    
    for (int i = 0; i < tf-1; ++i)
    {
        dp = (c_path[i+1] - c_path[i])*(Kd/(pow((Kd + c_path[i+1]),2)));
        dpdt_mean += (1./T)*dp *exp(((t_range[i]-t)/T));
    }

    return dpdt_mean;
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
    def("get_interpolation", get_interpolation);
    def("vorticity", vorticity);
    def("get_concentration", get_concentration);
    def("chemotaxis_all", chemotaxis_all);
    def("get_chemotaxis", get_chemotaxis);
    def("move_bacteria", move_bacteria);
    
};
