#!/usr/bin/env python

from math import *
from numpy import *
from pylab import *
import os.path
import matplotlib.patches as patches

params_input = sys.argv[1:]

fr_input = float(params_input[0])
taur_input = float(params_input[1])
alpha_input = float(params_input[2])

L = 2*pi
NP = 50000

filename = "FIG2.pdf"
fig = figure(figsize=(10., 4))


def get_data_pos(fr_input, taur_input, alpha_input, mtype):
    dir = "./data/data_FIG2/space-alpha%.3f-us%.3f-taur%.3f-%s" %(alpha_input, fr_input, taur_input, mtype)
    file_x = open(os.path.expanduser("%s" %dir))

    bacs = []
    for line in file_x:
        try:
            bacs.append(array([float(v) for v in line.split()]))
        except:
            pass
    bacs = array(bacs)
    
    return bacs



def get_data_vel(u, fr_input, taur_input, alpha_input, mtype):
    dir = "./data/data_FIG2/%s-alpha%.3f-us%.3f-taur%.3f-%s" %(u, alpha_input, fr_input, taur_input, mtype)
    file_x = open(os.path.expanduser("%s" %dir))

    u = []
    for line in file_x:
        try:
            u.append(array([float(v) for v in line.split()]))
        except:
            pass
    u = array(u)
    return u



for mtype in ["NN", "RR", "RT"]:
    if (mtype == "RR"):
        phi_use = pi
        rw_use = True
        aux = 3

    if (mtype == "RT"):
        phi_use = 1.23
        rw_use = True
        aux = 2

    if (mtype == "NN"):
        phi_use = 0
        rw_use = False
        aux = 1


    bacs = get_data_pos(fr_input, taur_input, alpha_input, mtype)
    omega = get_data_vel("vort", fr_input, taur_input, alpha_input, mtype)
    print("omega", omega)
    ax4 = subplot(1, 3, aux)
    print(omega.shape)
    xspace = linspace(0, L, omega.shape[0])
    yspace = linspace(0, L, omega.shape[0])
    X, Y = meshgrid(xspace, yspace)

    
    
    
    imshow(omega.T, origin = "lower", extent=[0,L,0,L], alpha = 0.75)

    for j in range(0, int(len(bacs[:, 0])/2)):
        patch = patches.Ellipse([bacs[j, 0], bacs[j, 1]], width = .145,  height = 0.0045, angle = (bacs[j, 2]/pi)*180, color = "black")
        ax4.add_patch(patch)

    if (aux > 1):
        ax4.yaxis.set_ticks([])
        ax4.set_ylabel(r" ")
    else:
        ax4.set_ylabel(r"$y$", fontsize=18)
    ax4.xaxis.set_ticks([0, L], [r"$0$", r"L"], fontsize=18)
    ax4.yaxis.set_ticks([0, L], [r"$0$", r"L"], fontsize=18)

    ax4.set_xlabel(r"$x$", fontsize=18)
    

    
ax4.annotate(r"(a)",  xycoords='figure fraction', xy= (0.05, 0.93), fontsize=18)
ax4.annotate(r"(b)",  xycoords='figure fraction', xy= (0.35, 0.93), fontsize=18)
ax4.annotate(r"(c)",  xycoords='figure fraction', xy= (0.67, 0.93), fontsize=18)

ax4.annotate(r"Simple swimmers",  xycoords='figure fraction', xy= (0.1, 0.82), fontsize=14)
ax4.annotate(r"Run-and-tumbling swimmers",  xycoords='figure fraction', xy= (0.39, 0.82), fontsize=14)
ax4.annotate(r"Run-reversing swimmers",  xycoords='figure fraction', xy= (0.73, 0.82), fontsize=14)



subplots_adjust(bottom=0.2, right=0.98, left=0.07, top = 0.8, wspace = 0.15, hspace = 0.9)
plt.savefig(filename)
