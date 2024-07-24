from ase.io import read, write
import sys
import os
import numpy as np

# 删除之前的输出文件（如果存在）
if os.path.isfile("ase_out.xyz"):
    os.remove("ase_out.xyz")

# 读取轨迹文件
traj = read(sys.argv[1], index=":")

# 初始化变量
max_avg_force_atoms = None
min_avg_force_atoms = None
max_avg_force = -np.inf
min_avg_force = np.inf

# 遍历所有结构，计算平均受力，找出最大和最小平均受力的结构
for atoms in traj:
    # 检查原子数量和元素种类
    if len(atoms) > 8 and all(element in atoms.get_chemical_symbols() for element in ['Ge', 'Sb', 'Te']):
        forces = atoms.get_forces()
        avg_force = np.mean(np.linalg.norm(forces, axis=1))
        
        if avg_force > max_avg_force:
            max_avg_force = avg_force
            max_avg_force_atoms = atoms.copy()
        
        if avg_force < min_avg_force:
            min_avg_force = avg_force
            min_avg_force_atoms = atoms.copy()

# 将结果写入文件
if max_avg_force_atoms:
    write("ase_out.xyz", max_avg_force_atoms, append=False)

if min_avg_force_atoms:
    write("ase_out.xyz", min_avg_force_atoms, append=True)
