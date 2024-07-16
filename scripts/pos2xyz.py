# Usage: python pos2xyz.py POSCAR
import sys
from ase.io import read, write

structure = read(sys.argv[1])
supercell = structure * (1, 1, 1)
write("model.xyz", supercell)
print("Over!!!")
