#!/usr/bin/env python

"""
用法: python filter_forces_in_trajectory.py [trajectory_file]

功能描述:
该脚本从给定的轨迹文件中筛选出原子受力在指定范围内的结构，并将这些结构写入到 `ase_out.xyz` 文件中。

参数说明:
- [trajectory_file]: 输入的轨迹文件路径（可以是任何ASE支持的格式，如XYZ、VASP等）。

结果:
- `ase_out.xyz` 文件，包含筛选后的结构。
"""

from ase.io import read, write
import sys
import os

def main():
    # 检查命令行参数
    if len(sys.argv) != 2:
        print("用法: python filter_forces_in_trajectory.py [trajectory_file]")
        sys.exit(1)

    trajectory_file = sys.argv[1]

    # 删除之前的输出文件（如果存在）
    if os.path.isfile("ase_out.xyz"):
        os.remove("ase_out.xyz")

    # 读取轨迹文件
    traj = read(trajectory_file, index=":")

    # 遍历所有结构，筛选符合条件的结构
    for atoms in traj:
        forces = atoms.get_forces()
        # 检查受力范围
        if forces.max() < 15 and forces.min() > -15:
            write("ase_out.xyz", atoms, append=True)

    print("筛选完成，结果已保存至 'ase_out.xyz' 文件中。")

if __name__ == "__main__":
    main()
