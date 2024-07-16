# Usage: python cif2vasp.py mp-2612.cif
import sys
import os
from ase.io import read, write
from ase.build import sort
path = sys.argv[1]
name = os.path.splitext(path)[0]
structure = read(path)
write('./{}.vasp'.format(name), sort(structure), direct=True)
