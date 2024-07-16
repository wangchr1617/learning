from pynep.calculate import NEP
from ase.io.vasp import read_vasp
import numpy as np
import sys

calc = NEP("nep.txt")
data = read_vasp(sys.argv[1])
des = calc.get_property('descriptor', data)
print(des.shape)
np.savetxt('des.txt', des, delimiter=',', newline='\n', fmt='%.18f')
