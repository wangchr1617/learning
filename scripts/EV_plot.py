# Usage: python EV_plot.py ev.txt

import matplotlib.pyplot as plt
import numpy as np
import sys
data = np.loadtxt(sys.argv[1], delimiter=',', dtype=float)
vol = data[:,0]
dft = data[:,1]

for i in range(len(dft)):
    if dft[i] == min(dft):
        x = vol[i]
        y = dft[i]
string = "({:.2f}, {:.2f})".format(x,y)

plt.figure(figsize = (5,3))
plt.plot(np.array(vol), np.array(dft), linestyle='--', alpha=0.7, label="DFT")
plt.scatter(np.array(vol), np.array(dft), alpha=0.7)
plt.scatter(x, y, c='r')
plt.text(x, y+0.05, string, verticalalignment="baseline", horizontalalignment="center")
plt.xlabel(r'Volume / $Å^3$')
plt.ylabel('Energy per atom / eV')
plt.legend(loc="best")
plt.savefig('./ev.png', bbox_inches='tight')
