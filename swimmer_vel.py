#!/usr/bin/env python

from math import *
from numpy import *
from numpy.random import normal
from numpy.fft import *
from pylab import *
from modules_source.randomflow_mod import *
from modules_source.swimmers_mod import *
from modules_source.foodgrad_mod import *
from modules_source.chemotaxis_mod import *
#from matplotlib.colors import LogNorm

import numpy as np, math, matplotlib.patches as patches
from mpl_toolkits.axes_grid1 import ImageGrid


params_input = sys.argv[1:]
vels = float(params_input[0])

a_file = open("Data_swimmers2020/motility/dtime_vels%.2f_us1_rr.txt" %(vels), 'w')

dt = 0.01
time = int(200./dt)
ug = 0.
dc_u = 1.

grid = Grid(N = 150, L = 2*pi)
randflow = Random_flow(grid, tau = 2*pi)
solver = ode(move_bac).set_integrator('dopri5')


swimmers_rT = Swimmer(uf = randflow.u, fr = vels, alpha = 0.98, phi0 = 1.23)
bac_rT = Ensemble(grid, N_bac = 500, chemo = True)

swimmers_rR = Swimmer(uf = randflow.u, fr = vels, alpha = 0.98, phi0 = pi)
bac_rR = Ensemble(grid, N_bac = 500 , chemo = True)

chemo = Chemotaxis(alpha=600, Kd = 1., sigma_rot = 0.01)
agg = Source(grid)
C = init_C(grid, agg, dc = dc_u)

for i in range(time):
   if (i == 0):  # only gets the flow (we need u(t0) and u(t1) for interpolation)
      ux, uy, eta = randflow.get_vel(grid.delta)
      continue

   ux_p, uy_p, eta = randflow.get_vel(grid.delta)
   C = evolve_C(grid, C, ux_p, uy_p, ug, i*dt, dt)
  
   agg.integrate(ux, uy, ux_p, uy_p,  grid, dt, i)
   C = agg.influx_source(grid, C, dt, dc_u, 0.01)
   randflow.flow_evolve(dt, grid)
   ux, uy = ux_p, uy_p

for i in range(time):

   ux_p, uy_p, eta = randflow.get_vel(grid.delta)

   C = evolve_C(grid, C, ux_p, uy_p, ug, i*dt, dt)

   for j in range(0, bac_rT.N_bac):
        chemotaxis(bac_rT, C, grid, i*dt, dt, chemo, swimmers_rT, j)
   for j in range(0, bac_rR.N_bac):
        chemotaxis(bac_rR, C, grid, i*dt, dt, chemo, swimmers_rR, j)
            
   swimmer_integrate(dt, i, bac_rT, ux, uy, ux_p, uy_p, ug, grid, swimmers_rT, solver)
   swimmer_integrate(dt, i, bac_rR, ux, uy, ux_p, uy_p, ug, grid, swimmers_rR, solver)

   agg.integrate(ux, uy, ux_p, uy_p,  grid, dt, i)
   C = agg.influx_source(grid, C, dt, dc_u, 0.01)
   randflow.flow_evolve(dt, grid)
   ux, uy = ux_p, uy_p

   dist_bacRT = []
   dist_bacRR = []
   for j in range(0, bac_rT.N_bac):
      dist_j = []
      for k in range(0, agg.N_a):
         dist_j.append(sqrt((bac_rT.x[j] - agg.x_s[k])**2 + (bac_rT.y[j] - agg.y_s[k])**2))
      dist_bacRT.append(min(dist_j))
      
   for j in range(0, bac_rR.N_bac):
      dist_j = []
      for k in range(0, agg.N_a):
         dist_j.append(sqrt((bac_rR.x[j] - agg.x_s[k])**2 + (bac_rR.y[j] - agg.y_s[k])**2))
      dist_bacRR.append(min(dist_j))


   a_file.write("%f, %f, %f \n" %(i*dt, mean(dist_bacRT)/grid.L, mean(dist_bacRR)/grid.L))
      
      
a_file.close()
   
