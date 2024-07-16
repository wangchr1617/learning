import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import sys

from ase.calculators.mixing import SumCalculator
from ase.io import read
from ase.units import GPa
from calorine.calculators import CPUNEP
from calorine.tools import relax_structure
from dftd3.ase import DFTD3

plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'

def NEP_calculator(dftd3=False):
    if os.path.exists("./nep.txt"):
        if dftd3 == True:
            calculator = SumCalculator([CPUNEP("./nep.txt"), DFTD3(method="pbe", damping="d3bj")])
        else:
            calculator = CPUNEP("./nep.txt")
        return calculator
    else:
        print("The file nep.txt does not exist.")
        sys.exit(1)

structure = read("./POSCAR")
structure.calc = NEP_calculator()

data = []
for volsc in np.arange(0.9, 1.1, 0.01):
    s = structure.copy()
    s.cell *= volsc
    s.calc = NEP_calculator()
    relax_structure(s, constant_volume=True)
    data.append(dict(volume=s.get_volume(),
                     energy=s.get_potential_energy() / len(structure),
                     pressure=-np.sum(s.get_stress()[:3]) / 3 / GPa))
df = pd.DataFrame(data)

plt.figure(figsize=(5, 3))
plt.plot(df.volume, df.energy, 'o-', alpha=0.7, label='E-V')

plt.xlabel('Volume / $Å^3$')
plt.ylabel('Energy per atom / eV')
plt.legend(loc="best")
plt.savefig('./ev_calorine.png', bbox_inches='tight')
