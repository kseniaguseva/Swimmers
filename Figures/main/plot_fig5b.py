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


startcolor =  '#a9c4f5'  #a9c4f5' 
midcolor =    '#fdf7ee' 
endcolor =    '#e5334b'

own_cmap3 = matplotlib.colors.LinearSegmentedColormap.from_list( 'own2', [startcolor, midcolor, endcolor] )


dt = 0.01
time = int(50./dt)
ug = 0.
dc_u = 1.

grid = Grid(N = 150, L = 2*pi)
randflow = Random_flow(grid, tau = 2*pi)
solver = ode(move_bac).set_integrator('dopri5')


swimmers_rT = Swimmer(uf = randflow.u, fr = 1., alpha = 0.98, phi0 = 1.23)
bac_rT = Ensemble(grid, N_bac = 500, chemo = True)

swimmers_rR = Swimmer(uf = randflow.u, fr = 1., alpha = 0.98, phi0 = pi)
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
   # if (i == 0):  # only gets the flow (we need u(t0) and u(t1) for interpolation)
   #    ux, uy = randflow.get_vel(grid.delta)
   #    continue

   #ux_p, uy_p, eta = randflow.get_vel(grid.delta)

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
   print(i*dt, bac_rT.N_bac, bac_rR.N_bac)
        
   if (i*dt > 40.):
        if(i%50 == 0):
          filename = "./Figures/swim_%.3f" %(i*dt)
          filename = filename + ".pdf"
          fig = figure(figsize=(6,4))
          grid_fig = ImageGrid(fig, 111,          # as in plt.subplot(111)
                               nrows_ncols=(1,2),
                               axes_pad=0.2,
                               share_all=False,
                               cbar_location="bottom",
                               cbar_mode="single",
                               cbar_size="2%",
                               cbar_pad=0.3,)
  
          ax1 = grid_fig[0]
          ax2 = grid_fig[1]
          omega, scos, ssin = get_vorticity(grid, ux, uy)
          omega = shift_frame(ug, i*dt, grid, omega)
          
          for j in range(0, bac_rT.N_bac):
               patch = patches.Ellipse([bac_rT.x[j]/grid.delta, bac_rT.y[j]/grid.delta], width = 5,  height = 0.75,angle = (bac_rT.theta[j]/pi)*180, color = "black")
               patch.update({'facecolor': 'gray','linewidth': 0.25 })
               
               ax1.add_patch(patch)
          for j in range(0, bac_rR.N_bac):
               patch = patches.Ellipse([bac_rR.x[j]/grid.delta, bac_rR.y[j]/grid.delta], width = 5,  height = 0.75,angle = (bac_rR.theta[j]/pi)*180, color = "black")
               patch.update({'facecolor': 'gray', 'linewidth': 0.25 })
               ax2.add_patch(patch)
          im = ax1.imshow(C.T, interpolation='nearest',  origin="lower", cmap= own_cmap3)
          im = ax2.imshow(C.T, interpolation='nearest',  origin="lower", cmap= own_cmap3)
          ax2.cax.colorbar(im)
          ax2.cax.toggle_label(True)

          ax1.set_title("Run-and-tumble")
          ax2.set_title("Run-reverse")
          ax1.contour(omega.T,  linewidths = 0.1, colors= 'blue')
          ax2.contour(omega.T,  linewidths = 0.1, colors= 'blue')

          ax1.set_xticks([0, grid.N])
          ax1.set_xticklabels([r"$0$", r"L"])
          ax1.set_yticks([0, grid.N])
          ax1.set_yticklabels([r"$0$", r"L"])
          ax1.xaxis.set_tick_params(labelsize=12)
          ax1.yaxis.set_tick_params(labelsize=12)

          ax2.set_xticks([0, grid.N])
          ax2.set_xticklabels([r"$0$", r"L"])
          ax2.xaxis.set_tick_params(labelsize=12)
          ax2.yaxis.set_tick_params(labelsize=12)
          ax1.set_xlim([0, grid.N])
          ax2.set_xlim([0, grid.N])
          ax1.set_ylim([0, grid.N])
          ax2.set_ylim([0, grid.N])
          ax1.set_xlabel("x")
          ax2.set_xlabel("x")
          ax1.set_ylabel("y")
          ax1.annotate(r"C",  xycoords='figure fraction', xy= (0.5, 0.03))
          subplots_adjust(bottom=0.1, right=0.98, left=0.1, top = 0.98, wspace = 0.45)
          plt.savefig(filename)
