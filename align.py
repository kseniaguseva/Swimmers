#!/usr/bin/env python

from math import *
from numpy import *
from pylab import *
from module_homogeneos.randomflow_mod import *
from module_homogeneos.swimmers_mod import *
import os.path
import time

##################################
#   INITIALIZATION
##################################

dt = 0.001                       # step size for integration
time = int(50./dt)              # simulation time  (dec was 60)
ug = 0.                         # no vertical motion

# input running time (time between switching direction)
params_input = sys.argv[1:]
fr_use = float(params_input[0])
tau_r = float(params_input[1])
alpha_use = float(params_input[2])
phi_type = str(params_input[3])
if phi_type == "NN":
    phi = 0      # no switching direction
if phi_type == "RT":
    phi = 1.23   # run and tumble
if phi_type == "RR":
    phi  = pi    #  run and reverse 

# set spatial grid and the flow field
grid = Grid(N = 250, L = pi)

randflow = Random_flow(grid, tau = 2*pi)

# sets the ode solver (use a function written in C)
solver = ode(move_bac_callC).set_integrator('dopri5')

##################################
# Parameters for bacterial motion
##################################
#  fr -- swimming speed u_s in relation to the velocity of the flow uf
#  alpha -- sets the shape of bacteria
#  phi -- sets the motility type
#  tau_r  -- running time, time period between direction switches
##################################

swimmers_rt = Swimmer(uf = randflow.u, fr = fr_use, alpha = alpha_use, phi0 = phi, tau_run = tau_r)
bac_rt = Ensemble(grid, N_bac = 10000)

##################################
#   TIME EVOLUTION
##################################

u_alignT = []
u_alignT_sqr = []
u_alignT_angle = []
for i in range(time):
    if (i == 0):  # only gets the flow (we need u(t0) and u(t1) for interpolation)
        ux, uy = randflow.get_vel(grid.delta)
        continue

    ux_p, uy_p = randflow.get_vel(grid.delta)  # for t > t0, gets u(t1)

    swimmer_integrate_periodic(dt, i, bac_rt, 0.1, ux, uy, ux_p, uy_p, ug, grid,
                                   swimmers_rt, solver, rw = True) # rw needs to be true for RR  and RT

    randflow.flow_evolve(dt, grid)
    ux, uy = ux_p, uy_p

    ######################################
    # computes <p.u>
    ######################################
    if(i*dt >= 2.):                  # leaves a transient out
        if(i%100 == 0):             # gets average at different time points
            u_align = []
            u_align_sqr = []
            u_align_angle = []
            for j in range(0, bac_rt.N_bac):
                ux_i = ux[int(bac_rt.x[j]/grid.delta), int(bac_rt.y[j]/grid.delta)]
                uy_i = uy[int(bac_rt.x[j]/grid.delta), int(bac_rt.y[j]/grid.delta)]
                #u_align.append(dot([ux_i/sqrt(ux_i**2+uy_i**2), uy_i/sqrt(ux_i**2+uy_i**2)], [cos(bac_rt.theta[j]), sin(bac_rt.theta[j])]))
                u_align_sqr.append(dot([ux_i/sqrt(ux_i**2+uy_i**2), uy_i/sqrt(ux_i**2+uy_i**2)], [cos(bac_rt.theta[j]), sin(bac_rt.theta[j])])**2)
                u_align.append(abs(dot([ux_i/sqrt(ux_i**2+uy_i**2), uy_i/sqrt(ux_i**2+uy_i**2)], [cos(bac_rt.theta[j]), sin(bac_rt.theta[j])])))
                u_align_angle.append(arccos(dot([ux_i/sqrt(ux_i**2+uy_i**2), uy_i/sqrt(ux_i**2+uy_i**2)], [cos(bac_rt.theta[j]), sin(bac_rt.theta[j])])))

            u_alignT_sqr.append(mean(u_align_sqr))
            u_alignT_angle.append(mean(u_align_angle))
            u_alignT.append(mean(u_align))
        # if(i*dt > 2.):                  # leaves a transient out
    #     if(i%10 == 0):             # gets average at different time points
    #         u_alignT.append(mean(abs(u_align)))

##################################
#   OUTPUT
##################################
    
filename = "data_final/align-us%.3f-tau%.3f_alpha%.3f_%s" %(fr_use, tau_r, alpha_use, phi_type)
print(filename)
if os.path.exists(filename):
    os.remove(filename)
temp_file = open(filename, "w")
for i in range(0, len(u_alignT)):
    print(u_alignT[i], u_alignT_sqr[i], u_alignT_angle[i], end="\n", file= temp_file)
   
temp_file.close()    
        
