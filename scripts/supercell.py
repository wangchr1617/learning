# Usage: python supercell.py mp-2612.vasp 2 2 2
import sys
import os
from ase.io import read, write
from ase.build import sort
path = sys.argv[1]
sx = eval(sys.argv[2])
sy = eval(sys.argv[3])
sz = eval(sys.argv[4])
name = os.path.splitext(path)[0]
structure = read(path) * (sx, sy, sz)
write('./{}_SC{}{}{}.vasp'.format(name,sx,sy,sz), sort(structure), direct=True)
