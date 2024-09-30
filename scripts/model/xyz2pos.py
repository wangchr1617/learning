# Usage: python xyz2pos.py sisb.xyz 8.6298163e+01 -5.9209473e-01  2.4067875e+00  4.3722048e+01  7.4572930e+01  9.5032241e-01  4.1447161e+01  2.4895231e+01  7.1588019e+01
from ase.io import read, write
from ase.io.vasp import write_vasp
import sys

def xyz_to_vasp(input_file, output_file, lattice_vectors=None):
    atoms = read(input_file, format='xyz')
    if lattice_vectors:
        atoms.set_cell(lattice_vectors)
        atoms.set_pbc(True)
    write_vasp(output_file, atoms, direct=True, sort=True)
input_file = sys.argv[1]
output_file = 'output.vasp'
lattice_vectors = [(sys.argv[2], sys.argv[3], sys.argv[4]), (sys.argv[5], sys.argv[6], sys.argv[7]), (sys.argv[8], sys.argv[9], sys.argv[10])]
xyz_to_vasp(input_file, output_file, lattice_vectors)
