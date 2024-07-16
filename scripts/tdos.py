# Usage: python tdos.py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

with open('./TDOS.dat','r') as f:
    data = []
    for line in f.readlines():
        l = []
        for item in line.replace("\n", "").split(" "):
            if item != "":
                l.append(item)
        data.append(l)
    df = pd.DataFrame(data)

x = pd.to_numeric(df[0][1:].reset_index(drop=True), errors='coerce')
y = pd.to_numeric(df[1][1:].reset_index(drop=True), errors='coerce')
plt.plot(x, y, color="r", linewidth=1, label="TDOS")
plt.xlabel('Energy(eV)'), plt.ylabel('DOS')
plt.xlim(x.min(), x.max()), plt.ylim(0, y.max()*1.1)
plt.legend(loc='best')
plt.savefig('./TDOS.png')
