#!/usr/bin/env python

from math import *
from numpy import *
from numpy.random import normal
from numpy.fft import *
import numpy.random
from numpy.random import *
from scipy import fftpack
from . import libfluid
from . randomflow_mod import *
from scipy.integrate import ode

from sys import *

class Swimmer():
     def __init__(self, uf = 1, fr = 1., alpha = 1., tau_run = 1., phi0 = pi):
         self.us = uf * fr
         self.alpha = alpha
         self.phi0 = phi0
         self.tau_run = tau_run

class Ensemble():
     def __init__(self, grid, N_bac = 1, chemo = False):
          self.N_bac = N_bac
          self.x = numpy.random.rand(self.N_bac)*grid.L
          self.y = numpy.random.rand(self.N_bac)*grid.L #zeros((self.N_bac))
          self.theta = numpy.random.rand(self.N_bac)*2*pi
          if(chemo == True):
               self.c_path = []
               self.c_last = zeros((self.N_bac, 2))
               self.c_total = []
               for i in range(0, self.N_bac):
                    self.c_path.append([])
               
     def boundary(self, grid):
          del_list = argwhere(self.y>grid.L).flatten()
          if (len(del_list) > 0):
               self.y = delete(self.y, del_list, 0)
               self.x = delete(self.x, del_list, 0)   
               self.theta = delete(self.theta, del_list, 0)
               self.c_last = delete(self.c_last, del_list, 0)
               for aux in del_list:
                    self.c_total.append(mean(self.c_path[aux]))
                    #print("deleted", self.c_total)
                    del self.c_path[aux]
                    
               self.N_bac = self.N_bac - len(del_list)
               
     def add_swimmers(self, grid, add_s = 10):
          x_new = numpy.random.rand(add_s)*grid.L
          theta_new = numpy.random.rand(add_s)*pi
          y_new = zeros((add_s))
          self.x = concatenate((self.x, x_new), axis=0)
          self.theta = concatenate((self.theta, theta_new), axis=0)
          self.y = concatenate((self.y, y_new), axis=0)
          for i in range(0, add_s):
               self.c_path.append([])
               add = zeros((1, 2))
               self.c_last = concatenate((self.c_last, add), axis = 0)
          self.N_bac = self.N_bac + add_s
     def add_cpath(self):
          #print(self.N_bac)
          for i in range(0, self.N_bac):
               self.c_path[i].append(self.c_last[i, 1])

               
def move_bac(t, x, ux, uy, ug, grid, swimmer, N_bac, tor):
    dtheta = []
    dx = []
    dy = []

    ux = shift_frame(ug, t, grid, ux)
    uy = shift_frame(ug, t, grid, uy)   
 
    omega, scos, ssin = get_vorticity(grid, ux, uy)
    for i in range(0, N_bac):
        x1 = x[i]
        x2 = x[N_bac + i]
        theta = x[2*N_bac + i]
            
        xn = grid.coordinates_grid(x1)
        yn = grid.coordinates_grid(x2)

        ux_pos = ux[xn, yn]
        uy_pos = uy[xn, yn]

        dtheta.append(omega[xn, yn] + swimmer.alpha*(scos[xn,yn]*cos(2*theta) - ssin[xn,yn]*sin(2*theta)))


        if(tor != None):
          tor.evolve_torque(omega[xn, yn], swimmer.alpha*(scos[xn,yn]*cos(2*theta) - ssin[xn,yn]*sin(2*theta)))
       
        dx.append(ux_pos + swimmer.us*cos(theta))
        dy.append(uy_pos + swimmer.us*sin(theta) + ug)
 
        
    return array([dx, dy, dtheta]).flatten()



def move_bac_callC(t, x, ux, uy, ux_p, uy_p, ug, grid, swimmer, N_bac, i, dt):
     dtheta = zeros(N_bac)
     dx = zeros(N_bac)
     dy = zeros(N_bac)

     ux = shift_frame(ug, t, grid, ux)
     uy = shift_frame(ug, t, grid, uy)   

     #print(t - (i-1)*dt)
     ux_u = ux + (t - (i-1)*dt)*(ux_p - ux)/dt
     uy_u = uy + (t - (i-1)*dt)*(uy_p - uy)/dt
     
     omega, scos, ssin = get_vorticity(grid, ux_u, uy_u)
     libfluid.move_bacteria(N_bac, ux_u, uy_u, omega, scos, ssin, x, dx, dy,
                            dtheta, grid.delta, grid.N, swimmer.alpha,
                            swimmer.us, ug)

     return array([dx, dy, dtheta]).flatten()


def swimmer_integrate(dt, i, bac, ux, uy, ug, grid, swimmer, tor = None):
     solver = ode(move_bac).set_integrator('dopri5')
     x0 = array([bac.x, bac.y, bac.theta]).flatten()
     solver.set_initial_value(x0, (i-1)*dt).set_f_params(ux, uy, ug, grid,
                                                         swimmer, bac.N_bac,
                                                         tor)
     solver.integrate(dt*i)

     bac.x = solver.y[:bac.N_bac]%grid.L
     bac.y = solver.y[bac.N_bac:2*bac.N_bac]
     bac.theta = solver.y[2*bac.N_bac:]%(2*pi)
     
     bac.boundary(grid)
     
def swimmer_integrate_periodic(dt, i, bac, eps, ux, uy, ux_p, uy_p, ug, grid,
                               swimmer, solver, rw = False):

     x0 = array([bac.x, bac.y, bac.theta]).flatten()
     solver.set_initial_value(x0, (i-1)*dt).set_f_params(ux, uy, ux_p, uy_p, ug,
                                                     grid, swimmer, bac.N_bac,
                                                     i, dt)
     solver.integrate(dt*i)

     bac.x = solver.y[:bac.N_bac]%grid.L
     bac.y = solver.y[bac.N_bac:2*bac.N_bac]%grid.L
     bac.theta = solver.y[2*bac.N_bac:]%(2*pi)
     if (rw == True):
          event = rand(bac.N_bac)
          for aux in range(bac.N_bac):
               if(event[aux] < dt/swimmer.tau_run):
                    angle = normal(0, eps)
                    bac.theta[aux] +=  swimmer.phi0 + angle
