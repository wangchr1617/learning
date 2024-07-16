import numpy as np
import os
import sys

from ase.calculators.mixing import SumCalculator
from ase.constraints import ExpCellFilter
from ase.io import read, write
from ase.optimize import BFGS
from ase.units import GPa
from calorine.calculators import CPUNEP
from dftd3.ase import DFTD3

def NEP_calculator(dftd3=False):
    if os.path.exists("./nep.txt"):
        if dftd3:
            calculator = SumCalculator([CPUNEP("./nep.txt"), DFTD3(method="pbe", damping="d3bj")])
        else:
            calculator = CPUNEP("./nep.txt")
        return calculator
    else:
        print("The file nep.txt does not exist.")
        sys.exit(1)

def relax(structure, pressure=0.0, maxstep=0.1, eps=None, max_step=None):
    structure.calc = NEP_calculator()
    mask = [True, True, True, True, True, True]
    ucf = ExpCellFilter(structure, scalar_pressure=pressure*GPa, mask=mask, constant_volume=False)
    gopt = BFGS(ucf, maxstep=maxstep)
    gopt.run(fmax=eps, steps=max_step)
    return structure
 
structure = read("./POSCAR") 
f_max = 0.001
Relaxed_atoms = relax(structure, eps=f_max, max_step=1000)

atom_force = Relaxed_atoms.get_forces()
U_atom = Relaxed_atoms.get_potential_energy()
if np.max(atom_force) <= f_max:
    print("  free  energy   TOTEN  =", U_atom, "eV")
    print("                 Voluntary context switches:")

write('CONTCAR', Relaxed_atoms, vasp5=True, direct=True)
