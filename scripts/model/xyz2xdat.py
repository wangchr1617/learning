from ase.io import read, write
import sys
atoms = read(sys.argv[1], index=':')
write('XDATCAR', atoms, format='vasp-xdatcar')

