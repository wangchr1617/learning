from ase.io import read, write
from pynep.calculate import NEP
from ase.optimize import BFGS
import numpy as np
from ase.constraints import ExpCellFilter
from ase.units import GPa

def relax(atoms, calc, pressure=0.0, maxstep=0.1, eps=None, max_step=None, dim=None):
    atoms.set_calculator(calc)
    mask = [True, True, True, True, True, True]
    if (dim == 2):
        ucf = ExpCellFilter(atoms, scalar_pressure=pressure*GPa, mask=mask, constant_volume=True)
    else:
        ucf = ExpCellFilter(atoms, scalar_pressure=pressure*GPa, mask=mask, constant_volume=False)

    gopt = BFGS(ucf, maxstep=maxstep)
    gopt.run(fmax=eps, steps=max_step)
    return atoms
 
model = read("./POSCAR") 
NEPcalc = NEP('./nep.txt')
f_max = 0.001
Relaxed_atoms = relax(atoms=model, calc=NEPcalc, eps=f_max, max_step=1000, dim=3)

atom_force = Relaxed_atoms.get_forces()
U_atom = Relaxed_atoms.get_potential_energy()
if (np.max(atom_force <= f_max )):
    print("  free  energy   TOTEN  =",U_atom,"eV")
    print("                 Voluntary context switches:")

write('CONTCAR', Relaxed_atoms, vasp5=True, direct=True)
