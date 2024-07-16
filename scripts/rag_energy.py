# Usage: python rag_energy.py
# grep Lattice model.xyz | awk '{print $21}' | cut -d '=' -f 2 > dft_energy.txt
# grep -B 1 'Lattice' model.xyz | grep -v 'Lattice' | grep -v -- '--' > atoms_num.txt

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
plt.rcParams['font.size'] = 14
plt.rcParams["axes.unicode_minus"] = False
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
import warnings
warnings.filterwarnings('ignore')

def qbc():
    ene = np.loadtxt('./dft_energy.txt', usecols=(0)).reshape((-1,1))
    natoms = np.loadtxt('./atoms_num.txt', usecols=(0)).reshape((-1,1))
    data = np.concatenate((ene, natoms), axis=1)
    df = pd.DataFrame(data, columns=["ene","natoms"])
    df = df[df["natoms"] == 62]
    df["ene"] = np.array(df["ene"]) / np.array(df["natoms"])
    return df

def count(a, x):
    xlist = []
    for i in x:
        xlist.append(np.sum(a == i))
    return np.array(xlist)

def plot(df, step=20):
    ini = np.round(np.array(df['ene'])[::20], decimals=3)
    rag = np.round(np.array(df['ene'])[::5], decimals=3)
    x1 = np.unique(sorted(ini))
    print(max(x1), min(x1))
    y1 = count(ini, x1)
    x2 = np.unique(sorted(rag))
    print(max(x2), min(x2))
    y2 = count(rag, x2)
    plt.figure(figsize=(5,3))
    plt.plot(x1, y1, color='black', label='Not relaxed')
    plt.plot(x2, y2, color='red', label='Relaxed')
    plt.xlim(-4.25, -4.18)
    plt.xlabel('Energy per atom (eV/atom)')
    plt.ylabel('Density of population')
    plt.savefig('./energy_population.png', bbox_inches='tight')

df = qbc()
plot(df)



