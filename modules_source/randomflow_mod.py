#!/usr/bin/env python

from math import *
from numpy import *
from numpy.random import normal
from numpy.fft import *
import numpy.random
from numpy.random import *
from scipy import fftpack
from . import libfluid
from numpy import linalg as LA

############################################################################

class Grid():
     def __init__(self, N = 125, L = 2*pi):
         self.N = N
         self.L = L
         self.delta = L/(N-1)
         self.x_i = linspace(0, L, N)
         self.y_i = linspace(0, L, N)
         x_grid, y_grid = meshgrid(self.x_i, self.y_i, indexing = "ij")
         self.x_grid = x_grid
         self.y_grid = y_grid
         
     def coordinates_grid(self, x):
          return int(x/self.delta)%(self.N)


def shift_frame(ug, t, grid, u):
     ############################################
     #   galilean transform: u(x-u_gt, t) + u_g #
     ############################################
    
     y0 = int(ug*t/grid.delta)%(grid.N)   # motion of the center of the vortex

     u = roll(u, y0, axis=1)              # shift the flow
 
     return u

###########################################################################

class Random_flow():
    def __init__(self, grid, tau = 1., l = 1., u = 1., seed = 15):
        self.tau = tau
        self.l = l
        self.u = u
        
        self.f_range = fftfreq(grid.N, 1./(grid.N))
        self.f_i = zeros(len(self.f_range) * 2)
        self.fmin = self.f_range.min()
        for i in range(len(self.f_range)):
            self.f_i[int(self.f_range[i] - self.fmin)] = i
        numpy.random.seed(seed)

        self.eta = zeros((grid.N, grid.N)) + zeros((grid.N, grid.N)) * 1j
        self.alpha = get_alpha(self.eta.shape, self.f_i, self.f_range, self.fmin)
        nu, mu = meshgrid(self.f_range, self.f_range, indexing = "ij")
        self.Q = self.u * sqrt(pi) * (self.l**2) * exp(-(((2*pi*l/grid.L) **2) * (mu**2 + nu**2))/4)

        self.eta = (grid.N*grid.delta)*self.Q * self.alpha 
        self.ieta = array(real(ifftn(self.eta)))/(grid.delta**2)
        
        
    def flow_evolve(self, dt, grid):
        alpha_n = get_alpha(self.eta.shape, self.f_i, self.f_range, self.fmin)
        self.alpha[:] = self.alpha * (1 - (dt/self.tau)) + sqrt(2*dt/self.tau)*alpha_n
        self.eta[:] = (grid.N*grid.delta) * self.Q * self.alpha
        self.ieta = array(real(ifftn(self.eta)))/(grid.delta**2)

    def get_vel(self, delta):
        vx = zeros(self.ieta.shape)
        vy = zeros(self.ieta.shape)
        
        libfluid.get_vel(self.ieta, vx, vy, delta)
        
        return vx, vy, self.ieta


def get_alpha(shape, f_i, f_range, fmin):
    alpha = zeros(shape) + zeros(shape) * 1j
    alpha = normal(0, sqrt(0.5), shape) + normal(0, sqrt(0.5), shape) * 1j
    alpha[0, 0] = alpha[int(shape[0]/2), 0] = alpha[0, int(shape[1]/2)] = alpha[int(shape[0]/2), int(shape[1]/2)] = normal(0, 1.)

    ralpha = alpha.real
    ialpha = alpha.imag

    libfluid.get_alpha(ralpha, ialpha, f_i, f_range, float(fmin))

    alpha = ralpha + ialpha * 1j

    return alpha
    



def get_vorticity(grid, ux, uy):
     omega = zeros((grid.N, grid.N))
     scos = zeros((grid.N, grid.N))
     ssin = zeros((grid.N, grid.N))
     libfluid.vorticity(omega, scos, ssin, ux, uy, grid.delta)
     return omega, scos, ssin
    

    
        

