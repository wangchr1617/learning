# Author: crwang
# This script can calculate Lattice paramters from BM state equation
# To use it: ./bm.py data
# data: ./energy_all.sh > data

import math
import numpy as np
import sys
# a = scaling coefficients, E = energy of structures
a, E = np.loadtxt(sys.argv[1], usecols=(0,1), delimiter=',', unpack=True)
x = (a * 4.18577)**(-2) # change 2.8664 to lattice para of 1.0/POSCAR
p = np.polyfit(x, E, 3)
c0 = p[3]
c1 = p[2]
c2 = p[1]
c3 = p[0]
x1 = (math.sqrt(4*c2**2 - 12*c1*c3) - 2*c2)/(6*c3)
para = 1/math.sqrt(x1)
print('The final lattice parameter is: %s ' %(para))
