from ase.calculators.mixing import SumCalculator
from ase.io import read
from calorine.calculators import CPUNEP
from calorine.tools import get_elastic_stiffness_tensor, relax_structure
from dftd3.ase import DFTD3

import sys
import os
import numpy as np

def NEP_calculator(dftd3=True):
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
# relax_structure(structure, fmax=0.0001)
cij = get_elastic_stiffness_tensor(structure, clamped=False, epsilon=0.005)
with np.printoptions(precision=1, suppress=True):
    print(cij)
