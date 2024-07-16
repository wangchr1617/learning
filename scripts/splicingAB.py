# Usage: python *.py A.vasp B.vasp
from ase.io import read, write
from ase.build import sort
import sys
import numpy as np

file_A = sys.argv[1] # './Si.vasp'
file_B = sys.argv[2] # './Sb.vasp'
structure_A = sort(read(file_A))
structure_B = sort(read(file_B))
cell_A = structure_A.cell
cell_B = structure_B.cell

spacing = 2.0
new_cell_A = np.array(cell_A)
new_cell_A[0:2] = np.maximum(cell_A[0:2], cell_B[0:2])
new_cell_B = np.array(cell_B)
new_cell_B[0:2] = np.maximum(cell_A[0:2], cell_B[0:2])
new_cell_C = np.zeros((3, 3))
new_cell_C[0:2] = np.maximum(cell_A[0:2], cell_B[0:2])
new_cell_C[2] = cell_A[2] + cell_B[2] + [0, 0, spacing]

structure_A.set_cell(new_cell_A, scale_atoms=False)
structure_B.set_cell(new_cell_B, scale_atoms=False)

zmax_A = structure_A.get_positions()[:, 2].max()
zmin_B = structure_B.get_positions()[:, 2].min()
z_translation = zmax_A - zmin_B + spacing

structure_B.translate((0, 0, z_translation))

structure_C = structure_A + structure_B
structure_C.set_cell(new_cell_C, scale_atoms=False)

write('AB.vasp', structure_C)
