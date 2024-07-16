# Usage: python *.py 001.vasp 001_REV.vasp
from ase import Atom
from ase.io import read, write
from ase.build import sort
import sys

filein = sys.argv[1]
fileout = sys.argv[2]
structure = sort(read(filein))

# find zmax, zmin
zmax = structure.get_positions()[:, 2].max()
zmin = structure.get_positions()[:, 2].min()

# find atoms B
def region(structure, X):
    AtomList = structure.get_chemical_symbols()
    begin = AtomList.index(X)
    for i in range(begin, len(AtomList)):
        if AtomList[i] != X:
            end = i
            break
        else:
            end = len(AtomList)
    return begin, end
begin, end = region(structure, 'B')

# insert atoms Li
for i in range(begin, begin+4):
    x, y, z = structure.get_positions()[i]
    z1 = zmax + 3
    z2 = zmin - 3
    structure = structure + Atom('Li',(x, y, z1))
    structure = structure + Atom('Li',(x, y, z2))
write(fileout, structure)
