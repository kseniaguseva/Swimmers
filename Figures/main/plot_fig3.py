#!/usr/bin/env python

from math import *
from numpy import *
from pylab import *
import os.path
import subprocess
from pandas import *

filename = "FIG3.pdf"
fig = figure(figsize=(9,3))

dt = 0.001
tau_r = 1.
alpha = 0.98
us = 0.5
us_range = arange(0.1, 4., 0.3)
tau_range = arange(0.1, 20., 0.5)
alpha_range = arange(0.1, 1., 0.1)

m = ["o", "s", "d"]
mc = ["black", "tomato", "lightskyblue"]
mlabel = ["s", "s+rr", "s+rt"]



def get_data_noise(tau_r, alpha, us, eps, align, s):
    df = read_csv("./data/data_FIG3/noise-us%.3f-tau%.3f_alpha%.3f_eps%.3f_%s" %(us, tau_r, alpha, eps, s),  delimiter = ' ', header=None)
    df.columns = ['abs', 'sqr', 'angle']
    df["time"] = arange(0, len(df["angle"])*dt, dt)
    return mean(df["abs"].values), std(df["abs"].values)


def get_data(tau_r, alpha, us, align, s):
    df = read_csv("./data/data_FIG3/align-us%.3f-tau%.3f_alpha%.3f_%s" %(us, tau_r, alpha, s),  delimiter = ' ', header=None)
    df.columns = ['abs', 'sqr', 'angle']
    df["time"] = arange(0, len(df["angle"])*dt, dt)
    return mean(df["abs"].values), std(df["abs"].values)


def get_data_angle(tau_r, alpha, us, align, s):
    df = read_csv("./data/data_FIG3/align-us%.3f-tau%.3f_alpha%.3f_%s" %(us, tau_r, alpha, s),  delimiter = ' ', header=None)
    df.columns = ['abs', 'sqr', 'angle']
    df["time"] = arange(0, len(df["angle"])*dt, dt)
    return mean(df["angle"].values*180/pi), std(df["angle"].values*180/pi)



def get_data_space(fr_input, taur_input, alpha_input, mtype):
    dir = "./data/data_FIG3/al-alpha%.3f-us%.3f-taur%.3f-%s" %(alpha_input, fr_input, taur_input, mtype)

    file_x = open(os.path.expanduser("%s" %dir))
    hist = []
    for x in file_x:
        try:

            for v in x.split()[0:]:
                hist.append(float(v))

        except:
            pass
    hist = array(hist)
    return hist


########################################################################################################
subplot(131)

NP = 50000
aux = 0
alpha_a = [0.3, 0.3, 0.3]
for mtype in ["NN", "RR", "RT"]:
    if (mtype == "RR"):
        phi_use = pi
        rw_use = True
    if (mtype == "RT"):
        phi_use = 1.23
        rw_use = True
    if (mtype == "NN"):
        phi_use = 0
        rw_use = False

    a = get_data_space(us, tau_r, alpha, mtype)
    hist, bins = histogram(a, bins = linspace(0, pi, 30))
    plot(bins[1:], hist/NP, "-",  color = mc[aux], label = mlabel[aux])
    bar(bins[1:], hist/NP, width = 0.1,  color = mc[aux], alpha = alpha_a[aux])

    aux += 1

annotate(r"(a)",  xycoords='figure fraction', xy= (0.01, 0.92), fontsize=15)
gca().set_xlabel(r"$\theta_u$", fontsize = 15)
gca().set_ylabel(r"$N/N_p$", fontsize = 15)
legend(loc = "best")

########################################################################################################
########################################################################################################
subplot(133)
c = ["black", "gray", "lightgray"]
aux = 0
for alpha_use in [0.98, 0.5, 0.25]:
    alignNN = []
    stdNN = []
    for us_use in us_range:
        m, s = get_data(tau_r, alpha_use, us_use, [], "NN")
        alignNN.append(m)
        stdNN.append(s)

    alignNN = array(alignNN)
    stdNN = array(stdNN)
    plot(us_range, alignNN, markeredgecolor = "black", marker = "o", markersize = 4, color = c[aux], label = r"$\alpha = %.2f$" %alpha_use)
    plt.fill_between(us_range, alignNN - stdNN, alignNN + stdNN, alpha = 0.2, color = c[aux])
    aux +=1

annotate(r"(b)",  xycoords='figure fraction', xy= (0.34, 0.92), fontsize=15)
gca().set_xlabel(r"$u_s/u_f$", fontsize = 15)
gca().set_ylabel(r"$\overline{\left<|\cos(\theta_u)|\right>}$", fontsize = 15)
legend(loc = "best")

########################################################################################################
########################################################################################################
subplot(132)
c = ["black", "gray", "lightgray"]
aux = 0
for alpha_use in [0.98, 0.5, 0.25]:
    alignNN = []
    stdNN = []
    for us_use in us_range:
        m, s = get_data_angle(tau_r, alpha_use, us_use, [], "NN")
        alignNN.append(m)
        stdNN.append(s)

    alignNN = array(alignNN)
    stdNN = array(stdNN)
    plot(us_range, alignNN, markeredgecolor = "black", marker = "o", markersize = 4, color = c[aux], label = r"$\alpha = %.2f$" %alpha_use)
    plt.fill_between(us_range, alignNN - stdNN, alignNN + stdNN, alpha = 0.2, color = c[aux])
    aux +=1


annotate(r"(c)",  xycoords='figure fraction', xy= (0.67, 0.92), fontsize=15)

gca().set_xlabel(r"$u_s/u_f$", fontsize = 15)
gca().set_ylabel(r"$\overline{\left<\theta_u\right>}$", fontsize = 15)


subplots_adjust(bottom=0.2, right=0.98, left=0.1, top = 0.9, hspace = 0.45, wspace = 0.45)
plt.savefig(filename)
