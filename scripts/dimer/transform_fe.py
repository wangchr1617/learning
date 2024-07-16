#!/usr/bin/env python
#Reset the first item of the energy data in
#fe.dat to zero and change the unit to eV

import numpy as np
import getopt
import sys

opts,args = getopt.getopt(sys.argv[1:],'-h',['help'])
for opt_name,opt_value in opts:
    if opt_name in ('-h','--help'):
        print('''#Reset the first item of the energy data in
#fe.dat to zero and change the unit to eV''')
        exit()

def read_from_fe(text):
    fp1 = open(text, mode='r')#,encoding='gbk')
    fp1.readline()
    data = fp1.readlines()
    A = [[] for i in range(len(data))]
    for i in range(len(data)):
        A[i][:] = list(map(float, data[i].strip('\n').split()))
    fp1.close
    return A

f = read_from_fe("fe.dat")
f1 = [[] for i in range(len(f))]
for i in range(len(f)):
    for n in range(2):
        f1[i].append(f[i][n])
#f[0].append(0)
f2 = np.array(f1)[:,0:2]
f2[:,1] = (f2[:,1] - f2[0,1]) * 27.2113838
np.savetxt("energy_eV.txt", np.array(f2))
