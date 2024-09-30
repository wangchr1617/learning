"""
用法: python filter_atoms.py [trajectory_file]

功能描述:
该脚本从给定的轨迹文件中筛选出符合条件的结构，
并将这些结构写入到 `ase_out.xyz` 文件中。

参数说明:
- [trajectory_file]: 输入的轨迹文件路径（可以是任何ASE支持的格式，如XYZ、VASP等）。

结果:
- `ase_out.xyz` 文件，包含筛选后的结构。
"""

from ase.io import read, write
import sys
import os

# 函数：过滤掉含有指定原子的结构
def filter_atoms(atoms, atom_symbol):
    """
    检查原子结构是否包含指定的元素符号。

    参数:
    atoms (ase.Atoms): 一个ASE的Atoms对象，表示原子结构。
    atom_symbol (str): 要过滤掉的原子符号（例如'Sb'）。

    返回值:
    bool: 如果结构中不含指定的原子，返回True；否则返回False。
    """
    return atom_symbol not in atoms.get_chemical_symbols()

# 函数：根据每原子的能量范围进行筛选
def filter_energies(atoms, energy_min, energy_max):
    """
    检查每原子的能量是否在指定范围内。

    参数:
    atoms (ase.Atoms): 一个ASE的Atoms对象，表示原子结构。
    energy_min (float): 每原子能量的最小允许值。
    energy_max (float): 每原子能量的最大允许值。

    返回值:
    bool: 如果每原子能量在指定范围内，返回True；否则返回False。
    """
    energy_per_atom = atoms.get_potential_energy() / len(atoms)
    return energy_min < energy_per_atom < energy_max
    
# 函数：根据力的范围进行筛选
def filter_forces(atoms, force_min, force_max):
    """
    检查原子结构的力是否在指定范围内。

    参数:
    atoms (ase.Atoms): 一个ASE的Atoms对象，表示原子结构。
    force_min (float): 原子受力的最小允许值。
    force_max (float): 原子受力的最大允许值。

    返回值:
    bool: 如果所有原子力在指定范围内，返回True；否则返回False。
    """
    forces = atoms.get_forces()
    return forces.max() < force_max and forces.min() > force_min

# 主函数：处理轨迹文件
def main():
    # 检查命令行参数
    if len(sys.argv) != 2:
        print("用法: python filter_atoms.py [trajectory_file]")
        sys.exit(1)

    # 获取输入的轨迹文件
    trajectory_file = sys.argv[1]

    # 删除之前的输出文件（如果存在）
    if os.path.isfile("ase_out.xyz"):
        os.remove("ase_out.xyz")

    # 读取轨迹文件
    traj = read(trajectory_file, index=":")

    # 遍历所有结构，筛选符合条件的结构
    for atoms in traj:
        # if filter_atoms(atoms, 'Sb') and filter_forces(atoms, 15.0) and filter_energies(atoms, -170.0, -150.0):
        if filter_energies(atoms, -170.0, -150.0):
            # 写入符合条件的结构到输出文件
            write("ase_out.xyz", atoms, append=True)

    print("筛选完成，结果已保存至 'ase_out.xyz' 文件中。")

if __name__ == "__main__":
    main()
