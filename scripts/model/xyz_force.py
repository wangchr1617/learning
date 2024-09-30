# Usage: python xxx.py 1.xyz 2.xyz

import matplotlib.pyplot as plt
import numpy as np
import sys
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'

def read_xyz(fileobj):
    lines = fileobj.readlines()
    data = []
    while len(lines) > 0:
        images = {}
        images["natoms"] = int(lines.pop(0))
        c0 = lines.pop(0) # Comment line
        c1 = c0.replace('forces', 'forces')
        c2 = c1.replace('Energy', 'energy')
        c3 = c2.replace('Virial', 'virial')
        c_list = c3.split(' ')
        # print(c_list) # 根据Properties位置决定下一行c_list[]的id, vasp_to_xyz.py提取出来的是9
        p = c_list[9].split('=')
        p_list = p[-1].split(':')
        d = {}
        pos = 0
        for i in range(int(len(p_list)/3)):
            dtype = p_list[3*i+1]
            width = p_list[3*i+2]
            begin = pos
            end = pos + eval(width)
            pos = end
            d[p_list[3*i]] = [dtype, begin, end, width]
        o = ['Properties=species',d["species"][0],d["species"][3],'pos',d["pos"][0],d["pos"][3],'forces',d["forces"][0],d["forces"][3]]
        comment = c3.replace(c_list[9], ':'.join(o))
        images["comment"] = comment
        symbols = []
        positions = []
        forces = []
        for _ in range(images["natoms"]):
            line = lines.pop(0) 
            symbol = line.split()[0]
            x, y, z = line.split()[d["pos"][1]:d["pos"][2]]
            fx, fy, fz = line.split()[d["forces"][1]:d["forces"][2]]
            symbol = symbol.lower().capitalize() # 元素符号首字母大写
            symbols.append(symbol)
            positions.append([float(x), float(y), float(z)])
            forces.append([float(fx), float(fy), float(fz)])
        images["symbols"] = symbols
        images["positions"] = positions
        images["forces"] = forces
        data.append(images)
    return data

def get_forcesarray(path):
    with open(path) as f:
        data = read_xyz(f)
        forces = []
        for images in data:
            for atoms in images["forces"]:
                for f in atoms:
                    forces.append(f)
    return np.array(forces)

path1 = sys.argv[1]
path2 = sys.argv[2]
f1 = get_forcesarray(path1)
f2 = get_forcesarray(path2)

plt.figure(figsize=(5,4))
plt.scatter(f1,f2)
plt.xlabel('{} force (eV/Å)'.format(path1))
plt.ylabel('{} force (eV/Å)'.format(path2))
plt.tight_layout()
plt.savefig('./force.png')

