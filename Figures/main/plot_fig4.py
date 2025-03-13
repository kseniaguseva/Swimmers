#!/usr/bin/env python

from math import *
from numpy import *
from pylab import *
import os.path
import subprocess
from pandas import *

filename = "FIG4.pdf"
fig = figure(figsize=(10,3))

dt = 0.001
tau_r = 1.
alpha = 0.98
us = 0.5
us_range = arange(0.1, 4., 0.3)
tau_range = arange(0.1, 12., 1.)
alpha_range = arange(0.1, 1., 0.1)

sm = ["d", "o", "s"]
mc = [ "lightskyblue", "black", "tomato"]
mlabel = ["s+rt", "s", "s+rr"]

def get_data(tau_r, alpha, us, align, s):
    df = read_csv("./data/data_FIG4/align-us%.3f-tau%.3f_alpha%.3f_%s" %(us, tau_r, alpha, s),  delimiter = ' ', header=None)
    df.columns = ['abs', 'sqr', 'angle']
    df["time"] = arange(0, len(df["angle"])*dt, dt)
    return mean(df["abs"].values), std(df["abs"].values)


########################################################################################################
subplot(131)

annotate(r"(a)",  xycoords='figure fraction', xy= (0.02, 0.92), fontsize=14)
annotate(r"(b)",  xycoords='figure fraction', xy= (0.37, 0.92), fontsize=14)
annotate(r"(c)",  xycoords='figure fraction', xy= (0.67, 0.92), fontsize=14)

aux = 0
for sv in ["RT", "NN", "RR"]:
    valign = []
    vstd = []
    for us_use in us_range:
        m, s = get_data(tau_r, alpha, us_use, [], sv)
        valign.append(m)
        vstd.append(s)

    valign = array(valign)
    vstd = array(vstd)

    plot(us_range, valign, marker = sm[aux], color = mc[aux],markeredgecolor = "black", markersize=4, label = mlabel[aux])
    plt.fill_between(us_range, valign - vstd, valign + vstd, color = mc[aux], alpha = 0.2)
    aux += 1

legend(loc = "best")
gca().set_xlabel(r"$u_s/u_f$", fontsize = 12)
gca().set_ylabel(r"$\overline{\left<|\cos(\theta_u)|\right>}$", fontsize = 12)

gca().set_xticks([0, 1, 2, 3, 4])

########################################################################################################
subplot(132)

aux = 0
for sv in ["RT", "NN", "RR"]:
    valign = []
    vstd = []

    for tau_r_use in tau_range:
        m, s = get_data(tau_r_use, alpha, us, [], sv)
        valign.append(m)
        vstd.append(s)
    valign = array(valign)
    vstd = array(vstd)


    plot(tau_range, valign, color = mc[aux], marker = sm[aux], markeredgecolor = "black", markersize=4, label = mlabel[aux])
    plt.fill_between(tau_range, valign - vstd, valign + vstd, color = mc[aux], alpha = 0.2)
    aux += 1


legend(loc = "best")

gca().set_xlabel(r"$\tau_0/t_f$", fontsize = 12)
gca().set_ylabel(r"$\overline{\left<|\cos(\theta_u)|\right>}$", fontsize = 12)


########################################################################################################
subplot(133)

aux = 0
for sv in ["RT", "NN", "RR"]:
    valign = []
    vstd = []
    for alpha in alpha_range:
        m, s = get_data(tau_r, alpha, us, [], sv)
        valign.append(m)
        vstd.append(s)
    valign = array(valign)
    vstd = array(vstd)


    plot(alpha_range, valign, color = mc[aux],  marker = sm[aux], markeredgecolor = "black", markersize=4, label = mlabel[aux])
    plt.fill_between(alpha_range, valign - vstd, valign + vstd, color = mc[aux], alpha = 0.2)
    aux += 1

gca().set_xlabel(r"$\alpha$", fontsize = 12)
gca().set_ylabel(r"$\overline{\left<|\cos(\theta_u)|\right>}$", fontsize = 12)

legend(loc = "best")

subplots_adjust(bottom=0.2, right=0.97, left=0.16, top = 0.95, hspace = 0.5, wspace = 0.65)
plt.savefig(filename)
