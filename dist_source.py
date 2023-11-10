
#!/usr/bin/env python

from math import *
from numpy import *
from numpy.random import normal
from numpy.fft import *
from pylab import *
import pandas as pd
import csv
from modules_source.randomflow_mod import *
from modules_source.swimmers_mod import *
from modules_source.foodgrad_mod import *
from modules_source.chemotaxis_mod import *


params_input = sys.argv[1:]
tau_input = float(params_input[0])
u_input = float(params_input[1])

f_save_round = open('./data/c_ensemble_tau%.2f_alpha0.00_us%.2f.csv' %(tau_input, u_input), 'w')

f_save_elip = open('./data/c_ensemble_tau%.2f_alpha0.98_us%.2f.csv' %(tau_input, u_input), 'w')

fields = ['mean', 'var', 'max', 'c_max', 'c_mean']
writer_round = csv.DictWriter(f_save_round, fieldnames=fields)
writer_elip = csv.DictWriter(f_save_elip, fieldnames=fields)
writer_round.writeheader()
writer_elip.writeheader()

dt = 0.01
time = int(500./dt)
ug = 0.0
init = 0 #int(0./dt)
dc_u = 0.1

grid = Grid(N = 250, L = 2*pi)
randflow = Random_flow(grid, tau = 2*pi)
c_alpha = [[], []]
aux_alpha = 0

agg = Source(grid)

solver = ode(move_bac).set_integrator('dopri5')

swimmers_r = Swimmer(uf = randflow.u, fr = tau_input, alpha = 0., phi0 = pi, tau_run = u_input)
bac_r = Ensemble(grid, N_bac = 500, chemo = True)

swimmers_e = Swimmer(uf = randflow.u, fr = tau_input, alpha = 0.98, phi0 = pi, tau_run = u_input)
bac_e = Ensemble(grid, N_bac = 500, chemo = True)

chemo = Chemotaxis(alpha=300, Kd = 1., sigma_rot = 0.1)

C, x_s, y_s = init_C(grid, dc = dc_u)
c_pos = []
for i in range(time):
   if (i == 0): # only gets the flow (we need u(t0) and u(t1) for
      #interpolation)
      ux, uy = randflow.get_vel(grid.delta)
      continue

   ux_p, uy_p = randflow.get_vel(grid.delta)

   C = evolve_C(grid, C, ux_p, uy_p, ug, i*dt, dt)
      
   for j in range(0, bac_r.N_bac):
      chemotaxis(bac_r, C, grid, i*dt, dt, chemo, swimmers_r, j)
   for j in range(0, bac_e.N_bac):
      chemotaxis(bac_e, C, grid, i*dt, dt, chemo, swimmers_e, j)

      
   swimmer_integrate(dt, i, bac_r, ux, uy, ux_p, uy_p, ug, grid,
                     swimmers_r, solver)
   swimmer_integrate(dt, i, bac_e, ux, uy, ux_p, uy_p, ug, grid,
                     swimmers_e, solver)

   
            
   agg.integrate(ux, uy, ux_p, uy_p,  grid, dt, i)
   C = agg.influx_source(grid, C, i*dt, dc_u, -0.001)
   randflow.flow_evolve(dt, grid)
   ux, uy = ux_p, uy_p

   c_aux_r = array(bac_r.c_path)
   c_aux_e = array(bac_e.c_path)
   #print(bac_rr.c_path[:][-1], len(c_aux[:, -1]), max(C.flatten()))
   #print("averages", mean(c_aux[:, -1]), mean(bac_rr.c_path[:][-1]))
   writer_round.writerow({'mean': mean(c_aux_r[:,-1]), 'var': var(c_aux_r[:,-1]), 'max': max(c_aux_r[:, -1]), 'c_max': max(C.flatten()), 'c_mean': mean(C.flatten())})
   writer_elip.writerow({'mean': mean(c_aux_e[:,-1]), 'var': var(c_aux_e[:,-1]), 'max': max(c_aux_e[:, -1]), 'c_max': max(C.flatten()), 'c_mean': mean(C.flatten())})
      

csvFile.close()
