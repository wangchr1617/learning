# Usage: python plt_msds.py MSD_1.dat MSD_2.dat
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'

plt.figure(figsize = (5,3), dpi=300)
max_x = []
max_y = []
for i in range(1, len(sys.argv)):
    data = np.loadtxt(sys.argv[i], skiprows=1, dtype=float)
    cols = ["Time", "x-MSD", "y-MSD", "z-MSD", "tot-MSD", "sqrt-MSD"]
    df = pd.DataFrame(data, columns=cols)
    label = sys.argv[i].split('.')[0]
    plt.plot(df[cols[0]], df[cols[4]], alpha=0.7, label=label)
    max_x.append(max(df[cols[0]]))
    max_y.append(max(df[cols[4]]))
plt.xlim(0, max(max_x))
plt.ylim(0, max(max_y)+1)
plt.xlabel("Time / fs")
plt.ylabel("MSD / Angstrom^2")
plt.title("Mean Square Displacement")
plt.legend()
plt.savefig("./msd.png", bbox_inches='tight')
