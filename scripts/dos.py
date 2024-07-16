# Usage: python dos.py Si Sb
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
import sys

element = sys.argv[1:]

with open('./TDOS.dat','r') as f:
    data = []
    for line in f.readlines():
        l = []
        for item in line.replace("\n", "").split(" "):
            if item != "":
                l.append(item)
        data.append(l)
    df0 = pd.DataFrame(data)
x0 = pd.to_numeric(df0[0][1:].reset_index(drop=True), errors='coerce')
y0 = pd.to_numeric(df0[1][1:].reset_index(drop=True), errors='coerce')
plt.plot(x0, y0, c="k", linewidth=1, label="TDOS")
plt.fill_between(x0, 0, y0, facecolor='grey', alpha=0.3)

def plt_pdos(path):
    label = os.path.splitext(path)[0]
    with open(path, 'r') as f:
        data = []
        for line in f.readlines():
            l = []
            for item in line.replace("\n", "").split(" "):
                if item != "":
                    l.append(item)
            data.append(l)
        df1 = pd.DataFrame(data)
    x1 = pd.to_numeric(df1[0][1:].reset_index(drop=True), errors='coerce')
    for i in range(4):
        y1 = pd.to_numeric(df1[i+1][1:].reset_index(drop=True), errors='coerce')
        if i == 0:
            y11 = y1
        else:
            y11 += y1
    plt.plot(x1, y11, alpha=0.7, linewidth=1, label=label)
for i in element:
    path = 'PDOS_{}.dat'.format(i)
    plt_pdos(path)

plt.xlabel('Energy(eV)'), plt.ylabel('DOS')
plt.xlim(x0.min(), x0.max())
plt.ylim(0, y0.max()*1.1)
plt.legend(loc='best')
plt.savefig('dos.png')
