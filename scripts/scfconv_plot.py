# Usage: python scfconv_plot.py 0

from ase.io import read
from ase.io.vasp import read_vasp, read_vasp_xdatcar
import numpy as np
import matplotlib.pyplot as plt
import sys
import warnings
warnings.filterwarnings('ignore')
plt.rcParams['font.size'] = 14
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'

idx = eval(sys.argv[1])
# Energy convergence
sc_steps_list = []
sc_energies_list = []
with open("OSZICAR", "r") as f:
    i = 0
    for line in f.readlines():
        if i > idx:
            break
        if "N " in line:
            sc_steps = []
            sc_energies = []
        elif "E0=" in line:
            sc_steps_list.append(sc_steps)
            sc_energies_list.append(sc_energies)
            i += 1
        else:
            parts = line.split()
            sc_steps.append(int(parts[1]))
            sc_energies.append(float(parts[2]))
plt.figure(figsize=(10, 6))
plt.grid(True)
plt.plot(sc_steps_list[idx], sc_energies_list[idx], '-o', label='Step {}'.format(idx))
plt.title('Energy of each self-consistency steps')
plt.xlabel('Self-consistency steps')
plt.ylabel('Energy(eV)')
plt.xlim(0, max(sc_steps_list[idx])+1)
plt.legend()
plt.savefig('scf.png', bbox_inches='tight')
