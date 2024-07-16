# Usage: python phonbydat.py band.dat

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'

with open(sys.argv[1]) as f:
    lines = f.readlines()
    comment = lines.pop(0)
    klabels = [eval(i) for i in lines.pop(0).strip().split(" ")[3:]]
    data = []
    band = []
    for line in lines:
        row = line.strip().split(" ")
        if "" in row:
            data.append(band)
            band = []
            continue
        band.append([eval(row[0]), eval(row[1])])
        
fig, ax = plt.subplots(figsize=(4.2, 3), dpi=140)
for i in range(len(data)):
    df = pd.DataFrame(data[i])
    for i in range(len(df)):
        x = df[:][0]
        y = df[:][1]
    ax.plot(x, y, c='r', lw=0.8)
ax.set_xlim(min(klabels), max(klabels))
ax.set_ylabel('Frequency (THz)')
ax.set_xticks(klabels)
# BAND_LABELS = [r'$\Gamma$', 'T', r'H$_2$', r'H$_0$', 'L', r'$\Gamma$', r'S$_0$', r'S$_2$', 'F', r'$\Gamma$']
# ax.set_xticklabels(BAND_LABELS)
for x in klabels:
    ax.axvline(x, c='0.5', ls="--", lw=0.8)
ax.axhline(y=0.0, c='0.5', ls="--", lw=0.8)
plt.tight_layout()
plt.savefig('./phon.png')