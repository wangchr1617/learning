# Usage: python exchangeABatoms.py POSCAR Ge Te 1 #001

from ase.io import read,write
from ase.build import sort
import os
import random
import sys

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

def rand_exchange(filename, A, B, n=1, idx=None):
    structure = sort(read(filename))
    name = os.path.splitext(filename)[0]
    Abegin, Aend = region(structure, A)
    Bbegin, Bend = region(structure, B)
    if Aend-Abegin < n or Bend-Bbegin < n:
        print("Error!!!")
        print("The atom number of element is less than the number you want delete!")
    
    if not os.path.exists('./exchange'):
        os.makedirs('./exchange')
    Alist = random.sample(range(Abegin, Aend), n)
    Blist = random.sample(range(Bbegin, Bend), n)
    if idx is None:
        idx = n
    for i in Alist:
        structure[i].symbol = B
    for i in Blist:
        structure[i].symbol = A
    write("./exchange/{}_{}.vasp".format(name, idx), sort(structure))
    #print("Sequence A:\n{}\nhas completed the exchange with \nSequence B:\n{}".format(Alist, Blist))

structure = sys.argv[1]
A = sys.argv[2]
B = sys.argv[3]
n = int(sys.argv[4])
idx = None
try:
    idx = sys.argv[5]
except:
    pass
rand_exchange(structure, A, B, n, idx)
