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


def advect_source(t, x,  ux, uy, ux_p, uy_p,  grid, dt, step, N_a):
     x_s = x[:N_a]
     y_s = x[N_a:]
     ux_u = ux + (t - (step-1)*dt)*(ux_p - ux)/dt
     uy_u = uy + (t - (step-1)*dt)*(uy_p - uy)/dt
     u_s = []
     v_s = []
     for i in range(0, N_a):
          pos_x = grid.coordinates_grid(x_s[i])
          pos_y = grid.coordinates_grid(y_s[i])
          u_s.append(ux_u[pos_x, pos_y])
          v_s.append(uy_u[pos_x, pos_y])
     dx_s = u_s 
     dy_s = v_s
     return array([dx_s, dy_s]).flatten()

class Source():
     def __init__(self, grid, N_a = 5):
          self.N_a = N_a
          self.x_s = numpy.random.rand(self.N_a)*grid.L
          self.y_s = numpy.random.rand(self.N_a)*grid.L
          self.solver =  ode(advect_source).set_integrator('dopri5')
          
                    
     def integrate(self, ux, uy, ux_p, uy_p,  grid, dt, i):
          x0 = array([self.x_s, self.y_s]).flatten()
         
          self.solver.set_initial_value(x0, dt*(i-1)).set_f_params(ux, uy, ux_p, uy_p,  grid,dt, i, self.N_a)
          self.solver.integrate(dt*i)
          self.x_s = self.solver.y[:self.N_a]
          self.y_s = self.solver.y[self.N_a:]
          
          self.x_s = self.x_s%grid.L
          self.y_s = self.y_s%grid.L
          
          
     def influx_source(self, grid, C, dt, dc_u,  rate = 0.01):
          C_new = C.copy()
          C_new -= C_new*rate*dt
          for i in range(0, self.N_a):
               pos_x = grid.coordinates_grid(self.x_s[i])
               pos_y = grid.coordinates_grid(self.y_s[i])
               C_new[pos_x, pos_y] = dc_u
          return C_new
     
     def dist_source(self, x, y):
          dist = []
          for i in range(0, self.N_a):
               dist.append(sqrt((self.x_s[i] - x)**2 + (self.y_s[i] - y)**2))
          return(min(dist))     
     

def init_C(grid, agg, dc = 1.):
     C = zeros((grid.N, grid.N))
     for i in range(0, agg.N_a):
          pos_x = grid.coordinates_grid(agg.x_s[i]) 
          pos_y = grid.coordinates_grid(agg.y_s[i])  
          C[pos_x, pos_y] = dc
         
     return C


def advect_back(ux, uy, grid, ug, t, dt):

     ux = shift_frame(ug, t, grid, ux)
     uy = shift_frame(ug, t, grid, uy)   

     x = grid.x_grid.copy() - ux*dt
     y = grid.y_grid.copy() - uy*dt - ug*dt
   
     return x, y

def evolve_C(grid, C, ux, uy, ug, t, dt):
     
     x, y = advect_back(ux, uy, grid, ug, t, dt)
     C_new = C.copy()
     
     libfluid.get_interpolation(grid.N, grid.N, x, y, 0., 0., grid.x_i,
                                grid.y_i, grid.delta, grid.delta, C, C_new)

     return C_new



