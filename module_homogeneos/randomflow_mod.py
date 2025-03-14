#!/usr/bin/env python

from math import *
from numpy import *
from numpy.random import normal
from numpy.fft import *
import numpy.random
import time
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
    def __init__(self, grid, tau = 1., l = 1., u = 1., seed = 10):
        self.tau = tau
        self.l = l
        self.u = u
        
        self.f_range = fftfreq(grid.N+1, 1./(grid.N+1))
        self.f_i = zeros(len(self.f_range) * 2)
        self.fmin = self.f_range.min()
        for i in range(len(self.f_range)):
            self.f_i[int(self.f_range[i] - self.fmin)] = i
            numpy.random.seed(int(time.time()))

        self.eta = zeros((grid.N+1, grid.N+1)) + zeros((grid.N+1, grid.N+1)) * 1j
        self.alpha = get_alpha(self.eta.shape, self.f_i, self.f_range, self.fmin)
        nu, mu = meshgrid(self.f_range, self.f_range)
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
        return vx, vy


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
    
# def get_strain_eign(vx, vy, grid, x, y):
#      i = int((x)/grid.delta)%grid.N
#      j = int((y)/grid.delta)%grid.N
#      ip = (i+1)%grid.N
#      im = (i-1)%grid.N
#      jp = (j+1)%grid.N
#      jm = (j-1)%grid.N

#      vxx = (vx[ip, j] - vx[im, j])/(2.*grid.delta)
#      vyy = (vy[i, jp] - vy[i, jm])/(2.*grid.delta)
#      vxy = (vx[i, jp] - vx[i, jm])/(2.*grid.delta)
#      vyx = (vy[ip, j] - vy[im, j])/(2.*grid.delta)
  
    
#      #E = array([[0.5*(vxx-vyy), 0.5*(vxy +vyx)],[0.5*(vyx +vxy), 0.5*(vyy-vxx)]])
#      E = array([[vxx, 0.5*(vxy +vyx)],[0.5*(vyx +vxy), vyy]])
     

#      T = E[0,0] + E[1, 1]
#      D = (E[0,0]*E[1, 1]) - (E[0,1]*E[1, 0])
#      L1 = (T/2.) + sqrt(((T**2)/4.) - D)
#      L2 = (T/2.) - sqrt(((T**2)/4.) - D)
#      v1 = zeros(2)
#      if(abs(E[0,1]) < 0.00001):  # E is diagonal
#           if(abs(E[0,0] - L1) > 0.0001):
#               v1 = [0, 1]
#           else:
#               v1 = [1, 0]  
#      else:
#           v1 = [1, -(E[0,0] - L1)/E[1,0]]

#      print(v1)
#      v1 = array(v1)
#      fc1 = 1./sqrt(v1[0]**2 + v1[1]**2)
#      v1 = v1*fc1
#      w, v = LA.eig(E)
    
#      return v1, w, v

def get_strain_eign(vx, vy, grid, x, y):
     i = int((x)/grid.delta)%grid.N
     j = int((y)/grid.delta)%grid.N
     ip = (i+1)%grid.N
     im = (i-1)%grid.N
     jp = (j+1)%grid.N
     jm = (j-1)%grid.N

     vxx = (vx[ip, j] - vx[im, j])/(2.*grid.delta)
     vyy = (vy[i, jp] - vy[i, jm])/(2.*grid.delta)
     vxy = (vx[i, jp] - vx[i, jm])/(2.*grid.delta)
     vyx = (vy[ip, j] - vy[im, j])/(2.*grid.delta)
  
    
     #E = array([[0.5*(vxx-vyy), 0.5*(vxy +vyx)],[0.5*(vyx +vxy), 0.5*(vyy-vxx)]])
     E = array([[vxx, 0.5*(vxy +vyx)],[0.5*(vyx +vxy), vyy]])
     
     # computing eigenvalues for 2x2 matrices
     
     T = E[0,0] + E[1, 1]
     D = (E[0,0]*E[1, 1]) - (E[0,1]*E[1, 0])
     L1 = (T/2.) + sqrt(((T**2)/4.) - D)
     L2 = (T/2.) - sqrt(((T**2)/4.) - D)
     v1 = zeros(2)
     if(abs(E[0,1]) < 0.00001):  # E is diagonal
          if(abs(E[0,0] - L1) > 0.0001):
              v1 = [0, 1]
          else:
              v1 = [1, 0]  
     else:
          v1 = [1, -(E[0,0] - L1)/E[1,0]]

     v1 = array(v1)
     fc1 = 1./sqrt(v1[0]**2 + v1[1]**2)
     v1 = v1*fc1

     # using LA library to compute the eigenvalues
     w, v = LA.eig(E)
     
     if(w[0] > w[1]):
         vp = v[:, 0]
         vs = v[:, 1]
     else:
         vp = v[:, 1]
         vs = v[:, 0]
    
     return v1, w, vp, vs


    

def alignment(ux, uy, bac, grid, ubins = 50):
     lambda1 = []
     lambda2 = []
     u_align = []

     for j in range(0, bac.N_bac):
          ux_i = ux[int(bac.x[j]/grid.delta), int(bac.y[j]/grid.delta)]
          uy_i = uy[int(bac.x[j]/grid.delta), int(bac.y[j]/grid.delta)]
          u_align.append(arccos(dot([ux_i/sqrt(ux_i**2+uy_i**2), uy_i/sqrt(ux_i**2+uy_i**2)], [cos(bac.theta[j]), sin(bac.theta[j])])))
    
          v1, w, vp, vs = get_strain_eign(ux, uy, grid, bac.x[j], bac.y[j])
          
          lambda1.append(arccos(dot(vp, [cos(bac.theta[j]), sin(bac.theta[j])])))
          
       
          lambda2.append(arccos(dot(vs, [cos(bac.theta[j]), sin(bac.theta[j])])))
          
       
          
     hist1, bins1 = histogram(lambda1, bins = ubins)
     hist2, bins2 = histogram(lambda2, bins = ubins)
     hist3, bins3 = histogram(u_align, bins = ubins)
     return hist1, bins1, hist2, bins2, hist3, bins3 

def alignment2(ux, uy, bac, grid, ubins = 50):
     lambda1 = []
     lambda2 = []
     u_align = []

     for j in range(0, bac.N_bac):
          ux_i = ux[int(bac.x[j]/grid.delta), int(bac.y[j]/grid.delta)]
          uy_i = uy[int(bac.x[j]/grid.delta), int(bac.y[j]/grid.delta)]
          u_align.append(arccos(dot([ux_i/sqrt(ux_i**2+uy_i**2), uy_i/sqrt(ux_i**2+uy_i**2)], [cos(bac.theta[j]), sin(bac.theta[j])])))

          v1, w, vp, vs = get_strain_eign(ux, uy, grid, bac.x[j], bac.y[j])

          lambda1.append(arccos(dot(vp, [cos(bac.theta[j]), sin(bac.theta[j])])))


          lambda2.append(arccos(dot(vs, [cos(bac.theta[j]), sin(bac.theta[j])])))



     return lambda1, lambda2, u_align
