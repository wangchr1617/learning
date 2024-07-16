# Usage: python rand_vacancy.py
from ase.io import read,write
from ase.build import sort
import os
import random

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

def rand_vacancy(filename, A, n):
    structure = sort(read(filename))
    begin, end = region(structure, A) 
    if end-begin < n:
        print("Error!!!")       
    Alist = random.sample(range(begin, end), n)
    print(Alist)
    del structure[Alist]
    write("./vac_{}.vasp".format(n), sort(structure))
    print("Over!!!")

rand_vacancy('./POSCAR', 'Ge', 150)
