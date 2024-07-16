# Usage: python EV_nep.py ev.txt nep.out

import matplotlib.pyplot as plt
import numpy as np
import sys
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'

path1 = sys.argv[1] 
path2 = sys.argv[2]
data = np.loadtxt(path1, delimiter=',', dtype=float)
nep = np.loadtxt(path2, delimiter=' ', dtype=float)
vol = data[:,0]
dft = data[:,1]

for i in range(len(dft)):
    if dft[i] == min(dft):
        x = vol[i]
        y = dft[i]
string = "({:.2f}, {:.2f})".format(x,y)

plt.figure(figsize = (5,3))
plt.plot(vol, nep, linestyle='-.', c='orange', alpha=0.7, label="NEP")
plt.plot(vol, dft, linestyle='--', alpha=0.7, label="DFT")
plt.scatter(np.array(vol), np.array(dft), alpha=0.7)
plt.scatter(x, y, c='red')
plt.text(x, y+0.1, string, verticalalignment="baseline", horizontalalignment="center")
plt.xlabel(r'Volume / $Å^3$')
plt.ylabel('Energy per atom / eV')
plt.legend(loc="upper right")
plt.savefig('./ev.png', bbox_inches = 'tight') 
