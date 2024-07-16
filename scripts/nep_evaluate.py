# Usage: python nep_evaluate.py nep_1.txt model.xyz 1
from ase.io import read
from pynep.calculate import NEP
import numpy as np
import sys
def evaluate(c, t, label):    
    calc = NEP(c)
    print(calc)
    traj = read(t, ':')
    e, f = [], []
    idx = 0
    for atoms in traj:
        num = len(atoms)
        atoms.set_calculator(calc)
        e.append([idx, atoms.get_potential_energy() / num])
        # f.append(atoms.get_forces().reshape(-1))
        idx += 1
    e = np.array(e)
    # f = np.concatenate(f)
    np.savetxt('nep_energy-{}.txt'.format(label), e, delimiter=' ', fmt="%.4f")
    # np.savetxt('nep_force-{}.txt'.format(label), f, delimiter=' ', fmt="%.4f")
evaluate(sys.argv[1], sys.argv[2], sys.argv[3])
