#!/usr/bin/env python

from math import *
from numpy import *
from pylab import *
import pandas as pd



filename = "FIG6.pdf"
                                       
fig = figure(figsize=(7, 2.5))

subplot(121)
shape_range = arange(0.1, 1, 0.1)
m_re = []
for shape in shape_range: 
    df = pd.read_csv("./data/data_FIG5/shape/dtime_sigma0,01_us1_shape%f.txt" %shape, sep=',', names=['time', 'RE'])
    
    re = df["RE"]
    t = df['time']
    
    m_re.append(mean(re[10000:]))
    
plot(shape_range, array(m_re)/0.22, "-d", label = r"rr", markeredgecolor = "black",color = "lightskyblue")

m_re = []
for shape in shape_range: 
    df = pd.read_csv("./data/data_FIG5/shape/dtime_sigma0,01_us1_shape%f_RT.txt" %shape, sep=',', names=['time', 'RT'])
    
    re = df["RT"]
    t = df['time']
    
    m_re.append(mean(re[10000:]))
    
plot(shape_range, array(m_re)/0.252, "-s", label = r"rt", markeredgecolor = "black", color =  "tomato")

gca().set_ylabel(r"$\frac{\overline{\left<d_{min}\right>}}{d^*}$", fontsize=16)
gca().set_xlabel(r"$\alpha$", fontsize=12)
legend(loc = "best")



subplot(122)
m_re = []
m_rr = []
vel_range = arange(0.1, 2.,0.1)
for vel in vel_range: 
    df = pd.read_csv("./data/data_FIG5/motility/dtime_vels%.2f_us1_rr.txt" %vel, sep=',', names=['time', 'RT', 'RR'])

    print(df['time'])
    re = df["RT"]
    rr = df["RR"]
    t = df['time']
    d = mean(re[0])
    print(d)  
    m_re.append(mean(re[10000:]))
    m_rr.append(mean(rr[10000:]))
  
plot(vel_range, array(m_rr)/d, "-d", label = r"rr", markeredgecolor = "black",color ="lightskyblue")
plot(vel_range, array(m_re)/d, "-s", label = r"rt", markeredgecolor = "black",color =  "tomato")

gca().set_ylabel(r"$\frac{\overline{\left<d_{min}\right>}}{d^*}$", fontsize=16)
gca().set_xlabel(r"$u_s/u_f$", fontsize=12)
legend(loc = "best")



annotate(r"(a)",  xycoords='figure fraction', xy= (0.01, 0.9), fontsize=12)
annotate(r"(b)",  xycoords='figure fraction', xy= (0.5, 0.9), fontsize=12)
subplots_adjust(bottom=0.2, right=0.98, left=0.15, top = 0.95, wspace = 0.5, hspace = 0.5)
plt.savefig(filename)
