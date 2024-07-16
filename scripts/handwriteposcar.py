# Usage: python handwriteposcar.py
from ase import Atoms
from ase.build import sort 
from ase.io import read, write
        
lattice_para = 6.015
angle = 89.954
atoms = Atoms('Ge4Te4', cell=[lattice_para, lattice_para, lattice_para, angle, angle, angle], pbc=(1, 1, 1), 
              scaled_positions=[(0.0,0.0,0.0), (0.0,0.5,0.5), (0.5,0.0,0.5), (0.5,0.5,0.0),
                                (0.5,0.0,0.0), (0.5,0.5,0.5), (0.0,0.0,0.5), (0.0,0.5,0.0)])
write("GeTe_cubic.vasp", sort(atoms), direct=True)
