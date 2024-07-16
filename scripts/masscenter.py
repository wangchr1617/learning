# -*- coding:utf-8 -*-

"""
@author: Jin-Cheng Liu
@file: masscenter
@time: 2/19/2017 10:03 PM
"""

import sys, math, linecache, time
import numpy as np


def atom_to_center_time():
    with open(fileName,'r') as f:
        lines = f.readlines()
        print('total lines of colored file:',len(lines))
    distance = []
    for i in range(0, loop_N):
        coord1_sum = [0, 0, 0]
        count_Au = 0
        for j in range(0, int(atomNum) + 2):
            tmpl = lines[i * (int(atomNum) + 2) + j].split(' ')
            while '' in tmpl:
                tmpl.remove('')
            if tmpl[0] == 'Au':
                count_Au += 1
                coord1_sum[0] += float(tmpl[1])
                coord1_sum[1] += float(tmpl[2])
                coord1_sum[2] += float(tmpl[3])

        coord_mass_center1 = list(map(lambda x: x/20, coord1_sum))
        c1 = np.array(tuple(coord_mass_center1))

        for j in range(0, int(atomNum) + 2):
            tmpl = lines[i * (int(atomNum) + 2) + j].split(' ')
            while '' in tmpl:
                tmpl.remove('')
            if tmpl[0] == 'Au':
                c2 = np.array([float(tmpl[1]),float(tmpl[2]),float(tmpl[3])])
                distance.append(str(np.linalg.norm(c1 - c2)))

    # print(distance)
    for i in distance:
        print(i)
    list2str = '\n'.join(distance)
    open(outfile, 'w').write(list2str)
    print('The distance between 2 cluster\'s mass center\
     are written to distance.txt')


print('#######  START  #######')
#fileName = 'newpos-CeO2.xyz'       # input filename
fileName = sys.argv[1]

atomNum = (linecache.getline(fileName,1)).strip()
print('total atom number:', atomNum)

#count the iteration number
with open(fileName,'r') as f:
    lines = f.readlines()
    print('total lines:',len(lines))

loop_N = ( len(lines) // (int(atomNum)+2) )
print('iteration number:',loop_N)

outfile = "distance.xyz"
atom_to_center_time()