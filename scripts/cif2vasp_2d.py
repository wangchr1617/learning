# Usage: python cif2vasp.py *.cif
import numpy as np
import sys
import os
from ase.io import read, write
from ase.build import sort

path = sys.argv[1]
name = os.path.splitext(path)[0]
structure = read(path)
# find cell
x, y, z = structure.get_cell()
zmax = structure.get_positions()[:, 2].max()
zmin = structure.get_positions()[:, 2].min()
z_3 = zmax-zmin
structure.set_cell([x, y, (0.0, 0.0, z_3)])
# find positions
pos = structure.get_positions()
new_pos = []
for i in range(len(pos)):
  x_i = pos[i][0]
  y_i = pos[i][1]
  z_i = pos[i][2] - zmin
  new_pos.append([x_i, y_i, z_i])
new_pos = np.array(new_pos)
structure.set_positions(new_pos)
# output
write('./{}.vasp'.format(name), structure, direct=False, sort=True)
