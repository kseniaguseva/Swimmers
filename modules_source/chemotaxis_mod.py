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


class Chemotaxis():
     def __init__(self, Kd = 30., T = 1, alpha = 660, sigma_rot = 0.25):
          self.Kd = Kd
          self.T = T
          self.alpha = alpha
          self.sigma_rot = sigma_rot

          
def get_Cpos(j, C, grid, bac, save = True):
     xn = int(bac.x[j]/grid.delta)%(grid.N)
     yn = int(bac.y[j]/grid.delta)%(grid.N)
     c_pos = C[xn, yn]
     if(save == True):
          bac.c_path[j].append(c_pos)     
     return c_pos

def chemotaxis(bac, C, grid, t, dt, chemo, swimmer, j):
     c_pos = get_Cpos(j, C, grid, bac)
     
     Pr = 0
     
     if (len(bac.c_path[j]) > 2.):
          ch = single_step_chemo(bac.c_path[j], c_pos, t, dt, chemo, swimmer)
     else:
          ch = swimmer.tau_run
        
    
     event = rand(1)
     if(ch > 1e-8):   # compute probability
          Pr = dt/ch
     else:            # the new tau_run is too small, turn
          Pr = 1.
                
     if (t == 0):
          Pr = 0.
        
     if event[0] < Pr:
          a = numpy.random.choice([-1, 1], 1)
          b = numpy.random.randn(1)*chemo.sigma_rot
          bac.theta[j] = (bac.theta[j] +(swimmer.phi0*a) + b)%(2*pi)
     #return c_pos


def single_step_chemo(c_path, c_pos, t, dt, chemo, swimmer):
     dp = (c_path[-1] - c_path[-2])*(chemo.Kd/(pow((chemo.Kd +c_path[-2]),2)))
     ch = swimmer.tau_run*exp(chemo.alpha*dp)
     return ch


def chemotaxis_new(bac, C, grid, t, dt, chemo, swimmer):
     libfluid.get_concentration(bac.N_bac, bac.c_last, C, bac.x,
                                bac.y, grid.delta, grid.N)
     bac.add_cpath()
     libfluid.chemotaxis_all(bac.c_last, bac.theta, chemo.Kd, chemo.alpha,
                             bac.N_bac, chemo.sigma_rot, swimmer.tau_run, dt,
                             swimmer.phi0)


