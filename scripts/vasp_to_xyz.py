# Usage: python vasp_to_xyz.py

from ase.io import read, write
import numpy as np
import sys
import os

if(len(sys.argv) == 2):
    label = sys.argv[1]
else:
    label= 'low'
os.system("find . -name vasprun.xml > xmllist")
os.system("if [ -f 'screen_tmp' ]; then rm screen_tmp; fi")
os.system("if [ -f 'NEP-dataset.xyz' ]; then rm NEP-dataset.xyz; fi")
for line in open('xmllist'):
    xml = line.strip('\n')
    print(xml)
    try:
        b = read(xml, index=":")
    except:
        b = read(xml.replace("vasprun.xml", "OUTCAR"), index=":")
        print(xml.replace("vasprun.xml", "OUTCAR"))
    # check convergence for each ionic step
    os.system("grep -B 1 E0 " + xml.replace('vasprun.xml','OSZICAR') + " | grep -E 'DAV|RMM' | awk '{if($2>=120) print 0; else print 1}' > screen_tmp")
    screen = np.loadtxt("screen_tmp")
    try:
        len(screen)
    except:
        screen=[screen]
    for ind, i in enumerate(screen):
        if(i == 1):
            xx, yy, zz, yz, xz, xy = -b[ind].calc.results['stress'] * b[ind].get_volume()
            b[ind].info['virial'] =  np.array([(xx, xy, xz), (xy, yy, yz), (xz, yz, zz)])
            del b[ind].calc.results['stress']
            b[ind].pbc = True
            b[ind].info['config_type'] = label
            write("NEP-dataset.xyz", b[ind], append=True)
    os.system("rm screen_tmp")
os.system("rm xmllist")

