#!/usr/bin/env python

from math import *
from numpy import *
from pylab import *


class Torque():
    def __init__(self):
         self.vort = []
         self.strain = []
    def evolve_torque(self, v, s):
        self.vort.append(v)
        self.strain.append(s)
    def out_torque(self, swimmer):
        np.savetxt("test/vort_%.3f.txt" %(swimmer.us), array(self.vort))
        np.savetxt("test/strain_%.3f.txt" %(swimmer.us), array(self.strain))