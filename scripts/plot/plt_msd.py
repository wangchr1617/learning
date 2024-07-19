# Usage: python plt_msd.py
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'

data = np.loadtxt('./MSD.dat', skiprows=1, dtype=float)
cols = ["Time", "x-MSD", "y-MSD", "z-MSD", "tot-MSD", "sqrt-MSD"]
df = pd.DataFrame(data, columns=cols)

# Column "sqrt-MSD" data not used
plt.figure(figsize = (5,3))
colors = ["forestgreen", "orange", "orangered", "chocolate", "dodgerblue", "blueviolet"]
for i in range(4):
    plt.plot(df[cols[0]], df[cols[i+1]], c=colors[i+1], label=cols[i+1])
plt.xlim(0, max(df[cols[0]]))
plt.ylim(0, max(df[cols[4]])+1)
plt.xlabel("t / fs")
plt.ylabel("MSD / Angstrom^2")
plt.title("Mean Square Displacement")
plt.legend()
plt.savefig("./msd.png", bbox_inches='tight')
