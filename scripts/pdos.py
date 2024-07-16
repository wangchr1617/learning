# Usage: python pdos.py PDOS_Sb.dat
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import sys

path = sys.argv[1]
target = os.path.splitext(path)[0]+'.png'
with open(path, 'r') as f:
    data = []
    for line in f.readlines():
        l = []
        for item in line.replace("\n", "").split(" "):
            if item != "":
                l.append(item)
        data.append(l)
    df = pd.DataFrame(data)

x = pd.to_numeric(df[0][1:].reset_index(drop=True), errors='coerce')
y_max = []
for i in range(4):
    y = pd.to_numeric(df[i+1][1:].reset_index(drop=True), errors='coerce')
    plt.plot(x, y, alpha=1, linewidth=1, label=df[i+1][0])
    y_max.append(y.max())
plt.xlabel('Energy(eV)'), plt.ylabel('DOS(states/eV)')
plt.xlim(x.min(), x.max()), plt.ylim(0, max(y_max)*1.1)
plt.legend(loc='best')
plt.savefig(target)
