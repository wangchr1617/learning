"""
用法: python analyze_forces_in_trajectory.py [trajectory_file]

功能描述:
该脚本用于读取原子轨迹文件，计算每个结构的平均受力，并找出平均受力最大的结构和最小的结构。
结果将写入到 `ase_out.xyz` 文件中。

参数说明:
- [trajectory_file]: 输入的轨迹文件路径（可以是任何ASE支持的格式，如XYZ、VASP等）。

结果:
- `ase_out.xyz` 文件，包含最大和最小平均受力的结构。
"""

from ase.io import read, write
import sys
import os
import numpy as np

def main():
    # 检查命令行参数
    if len(sys.argv) != 2:
        print("用法: python analyze_forces_in_trajectory.py [trajectory_file]")
        sys.exit(1)

    trajectory_file = sys.argv[1]

    # 删除之前的输出文件（如果存在）
    if os.path.isfile("ase_out.xyz"):
        os.remove("ase_out.xyz")

    # 读取轨迹文件
    traj = read(trajectory_file, index=":")

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
    if max_avg_force_atoms is not None:
        write("ase_out.xyz", max_avg_force_atoms, append=False)
        print(f"已将最大平均受力结构写入到 'ase_out.xyz'，平均受力为: {max_avg_force:.3f}")

    if min_avg_force_atoms is not None:
        write("ase_out.xyz", min_avg_force_atoms, append=True)
        print(f"已将最小平均受力结构追加到 'ase_out.xyz'，平均受力为: {min_avg_force:.3f}")

if __name__ == "__main__":
    main()
