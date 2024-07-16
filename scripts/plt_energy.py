#!/bin/python
# Author: crwang
# This script can input data and draw a curve. The figure will be saved as fig_name
# To use it: ./plt.py data_name fig_name

import matplotlib.pyplot as plt
import numpy as np
import sys
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'

x,y = np.loadtxt(sys.argv[1], delimiter=',', usecols=(0,1), unpack=True)
plt.xlabel('ENCUT (eV)')
plt.ylabel('Static energy (eV/atom)')
plt.plot(x, y, 'rs-', linewidth=2.0, label='k-mesh = 7×7×7')
plt.legend(loc='best')
plt.savefig('./' + sys.argv[2], bbox_inches='tight')
