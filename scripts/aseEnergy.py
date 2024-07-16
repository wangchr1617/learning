# Usage: python aseEnergy.py
from ase.calculators.emt import EMT
from ase.calculators.lj import LennardJones
from ase.io import read

atoms = read('POSCAR')
# atoms.set_calculator(EMT())
atoms.set_calculator(LennardJones())
ene = atoms.get_potential_energy()
print(ene)

