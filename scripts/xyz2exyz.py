# Usage: python xyz2exyz.py train.xyz

import sys

def read_xyz(fileobj):
    lines = fileobj.readlines()
    data = []
    while len(lines) > 0:
        images = {}
        images["natoms"] = int(lines.pop(0))
        c0 = lines.pop(0) # Comment line
        c1 = c0.replace('dft_forces', 'forces')
        c2 = c1.replace('dft_energy', 'energy')
        c3 = c2.replace('virial_stress', 'virial')
        c_list = c3.split(' ')
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

def write_exyz(fileobj, data):
    for images in data:
        natoms = images["natoms"]
        comment = images["comment"].rstrip()
        if '\n' in comment:
            raise ValueError('Comment line should not have line breaks.')
        fileobj.write('%d\n%s\n' % (natoms, comment))
        for i in range(images["natoms"]):
            s = images["symbols"][i]
            x, y, z = images["positions"][i]
            fx, fy, fz = images["forces"][i]
            fileobj.write('{:<2} {:>16.8f} {:>16.8f} {:>16.8f} {:>16.8f} {:>16.8f} {:>16.8f}\n'.format(s, x, y, z, fx, fy, fz))

path = sys.argv[1]
with open(path) as f:
    data = read_xyz(f)
    print("Number of structures: ",len(data))
with open('./out.xyz', 'a+') as f:
    write_exyz(f, data)
