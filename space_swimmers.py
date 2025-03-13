#!/usr/bin/env python

from math import *
from numpy import *
from pylab import *
from module_homogeneo.randomflow_mod import *
from module_homogeneo.swimmers_mod import *
import matplotlib.patches as patches
import os.path


params_input = sys.argv[1:]

fr_input = float(params_input[0])
taur_input = float(params_input[1])
alpha_input = float(params_input[2])


dt = 0.01
time = int(10./dt)
ug = 0.

for mtype in ["NN", "RR", "RT"]:
    if (mtype == "RR"):
        phi_use = pi
        rw_use = True
        aux = 6
    if (mtype == "RT"):
        phi_use = 1.23
        rw_use = True
        aux = 5
    if (mtype == "NN"):
        phi_use = 0
        rw_use = False
        aux = 4
    
    grid = Grid(N = 250, L = 2*pi)
    randflow = Random_flow(grid, tau = 2*pi)

    solver = ode(move_bac_callC).set_integrator('dopri5')
    swimmers_rt = Swimmer(uf = randflow.u, fr = fr_input, alpha = alpha_input, phi0 = phi_use, tau_run = taur_input)
    bac_rt = Ensemble(grid, N_bac = 50000)


    for i in range(time):
        if (i == 0):  # only gets the flow (we need u(t0) and u(t1) for interpolation)
            ux, uy = randflow.get_vel(grid.delta)
            continue

        ux_p, uy_p = randflow.get_vel(grid.delta)  # for t > t0, gets u(t1)


        swimmer_integrate_periodic(dt, i, bac_rt, ux, uy, ux_p, uy_p, ug, grid,
                                   swimmers_rt, solver, rw = rw_use)

        randflow.flow_evolve(dt, grid)
        ux, uy = ux_p, uy_p
        print(i*dt, mtype)
        
    hist1, bins1, hist2, bins2, hist3, bins3 = alignment(ux, uy, bac_rt, grid)
    
    filename = "data/al-alpha%.3f-us%.3f-taur%.3f-%s" %(alpha_input, fr_input, taur_input, mtype)
    if os.path.exists(filename):
        os.remove(filename)
    temp_file = open(filename, "w")
    for i in range(0, len(bins1)-1):
        print(bins1[i+1], hist1[i], bins2[i+1], hist2[i], bins3[i+1], hist3[i], sep=' ', end="\n", file= temp_file)
      
    temp_file.close()

    filename = "data/space-alpha%.3f-us%.3f-taur%.3f-%s" %(alpha_input, fr_input, taur_input, mtype)
    if os.path.exists(filename):
        os.remove(filename)
    temp_file = open(filename, "w")
    for i in range(0, bac_rt.N_bac):
        print(bac_rt.x[i], bac_rt.y[i], bac_rt.theta[i], sep=' ', end="\n", file= temp_file)
       

    temp_file.close() 


