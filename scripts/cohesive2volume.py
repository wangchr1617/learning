import numpy as np
import re

def calculate_cell_volume(xyz_file):
    with open(xyz_file, 'r') as file:
        lines = file.readlines()
        lattice_line = lines[1]
        lattice_match = re.search(r'Lattice="([^"]+)"', lattice_line)
        if lattice_match:
            lattice_params = lattice_match.group(1)
            lattice_vectors = np.array([float(x) for x in lattice_params.split()])
            vec_a = lattice_vectors[0:3]
            vec_b = lattice_vectors[3:6]
            vec_c = lattice_vectors[6:9]
            a = np.array(vec_a)
            b = np.array(vec_b)
            c = np.array(vec_c)
            volume = np.abs(np.dot(a, np.cross(b, c)))
            return volume
        else:
            raise ValueError("Lattice parameter not found in the .xyz file.")

volume = calculate_cell_volume('model.xyz')
natoms = 2

with open('cohesive.out', 'r') as f, open('ev.txt', 'w') as fout:
    for line in f:
        # 分割每行为scale和energy
        split_line = line.split()
        scale = float(split_line[0])
        energy = float(split_line[1])

        # 计算新的体积和能量
        new_volume = (scale ** 3) * volume
        new_energy_per_atom = energy / natoms

        # 写入新的体积和能量至输出文件
        fout.write(f"{new_volume:.10e},{new_energy_per_atom:.10e}\n")

print("The ev.txt file has been created with updated values.")


