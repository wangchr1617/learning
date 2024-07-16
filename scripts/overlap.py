from ase import Atoms
from ase.io import read, write
from scipy.spatial import cKDTree
import numpy as np

filename = "poly.xyz"
atoms = read(filename)
positions = atoms.get_positions()

tree = cKDTree(positions)
threshold = 1.5
close_pairs = tree.query_pairs(threshold)

atoms_to_remove = set()
for i, j in close_pairs:
    atom_to_remove = min(i, j)
    atoms_to_remove.add(atom_to_remove)
new_atoms = Atoms([atom for idx, atom in enumerate(atoms) if idx not in atoms_to_remove])

output_filename = "modified_poly.xyz"
write(output_filename, new_atoms)

print(f"Modified file saved as {output_filename}")

