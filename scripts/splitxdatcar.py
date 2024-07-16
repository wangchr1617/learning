# from ase import Atom
from ase.io import read, write
# from ase.build import sort
import sys

begin = eval(sys.argv[1])
end = eval(sys.argv[2])
try:
  step = eval(sys.argv[3])
except:
  step = None
traj = read('./XDATCAR', index=slice(begin, end, step))
print('The length of traj is: ', len(traj))
write('./XDATCAR_REV', traj, format='vasp-xdatcar')
