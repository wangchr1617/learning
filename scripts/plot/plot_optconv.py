"""
文件名: plot_optconv.py
运行方法: python plot_optconv.py ./opt/ (POSCAR_ini)
功能描述: 从 VASP 输出文件中读取数据并生成能量收敛、力收敛和初始结构变化的图表。
"""

from ase.io import read
from ase.io.vasp import read_vasp, read_vasp_xdatcar
import numpy as np
import matplotlib.pyplot as plt
import sys
import warnings

# 忽略警告信息
warnings.filterwarnings('ignore')

# 设置全局画图参数
plt.rcParams['font.size'] = 14
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'

def write_force(outdata, atoms):
    """
    从 OUTCAR 文件中提取最大力。

    参数:
    outdata (list): OUTCAR 文件的行列表。
    atoms (ase.Atoms): ASE Atoms 对象，用于约束信息。

    返回:
    list: 各优化步骤的最大力列表。
    """
    fmax_list = []
    f = []
    for x, line in enumerate(outdata):
        if "TOTAL-FORCE" in line:
            f = []
            for y in range(x + 2, len(outdata)):
                # 检查是否达到力的末尾
                if outdata[y][1] == "-":
                    break
                # 提取力数据
                f.append(outdata[y].split()[3:6])
            try:
                # 获取原子的约束信息
                c = atoms._get_constraints()
                indices_fixed = c[0].index
                for i in indices_fixed:
                    f[i] = [0, 0, 0]  # 对约束原子力设置为零
            except:
                pass
            fmax = 0
            # 计算最大力
            for i in f:
                fval = (float(i[0])**2 + float(i[1])**2 + float(i[2])**2)**(1./2.)
                if fval > fmax:
                    fmax = fval
            fmax_list.append(fmax)
    return fmax_list

def cartesian_to_fractional(structure):
    """
    将笛卡尔坐标转换为分数坐标。

    参数:
    structure (ase.Atoms): ASE Atoms 对象，包含原子坐标和晶胞信息。

    返回:
    np.array: 分数坐标数组。
    """
    cell = np.matrix(structure.cell)
    X = (cell.T)**-1
    cart_coords = np.matrix(structure.get_positions())
    frac_coords = (X * cart_coords.T).T
    return np.array(frac_coords)

def plot_convergence(path, initial_structure_path=None):
    """
    从指定路径读取数据并绘制收敛图。

    参数:
    path (str): 包含 VASP 输出文件的路径。
    initial_structure_path (str): 初始结构的路径，可选。
    """
    # 能量收敛图
    ion_steps = []
    energies = []
    try:
        with open(f"{path}/OSZICAR", "r") as f:
            for line in f.readlines():
                if "E0" in line:
                    parts = line.split()
                    ion_steps.append(int(parts[0]))
                    energies.append(float(parts[4]))
    except FileNotFoundError:
        print(f"Error: 文件 '{path}/OSZICAR' 未找到。")
        sys.exit(1)

    plt.figure(figsize=(10, 6))
    plt.grid(True)
    plt.plot(ion_steps, energies, '-o', label='Energy')
    plt.title('Energy Convergence')
    plt.xlabel('Ion steps')
    plt.ylabel('Energy (eV)')
    plt.xlim(0, max(ion_steps) + 1)
    plt.legend()
    plt.savefig('conv_ene.png', bbox_inches='tight')

    # 力收敛图
    try:
        out = open(f"{path}/OUTCAR").readlines()
        atoms = read(f"{path}/CONTCAR")
    except FileNotFoundError:
        print(f"Error: 文件 '{path}/OUTCAR' 或 '{path}/CONTCAR' 未找到。")
        sys.exit(1)

    forces = write_force(out, atoms)
    plt.figure(figsize=(10, 6))
    plt.grid(True)
    plt.plot(range(1, len(forces) + 1), forces, "-o", label="Force")
    plt.title('Force Convergence')
    plt.xlabel('Optimization Step')
    plt.ylabel('Max Force (eV/Å)')
    plt.xlim(0, len(forces) + 1)
    plt.legend()
    plt.savefig('conv_force.png', bbox_inches='tight')

    # 与初始结构的距离变化
    if initial_structure_path is None:
        initial_structure = read_vasp(f"{path}/POSCAR")
    else:
        print("警告！POSCAR 不是您定义的初始结构。")
        initial_structure = read_vasp(initial_structure_path)

    try:
        structures = read_vasp_xdatcar(f"{path}/XDATCAR", index=slice(None))
    except FileNotFoundError:
        print(f"Error: 文件 '{path}/XDATCAR' 未找到。")
        sys.exit(1)

    distances = []
    for structure in structures:
        rmsd = np.sqrt(((cartesian_to_fractional(structure) - cartesian_to_fractional(initial_structure)) ** 2).sum(axis=1).mean())
        distances.append(rmsd)

    plt.figure(figsize=(10, 6))
    plt.grid(True)
    plt.plot(range(1, len(distances) + 1), distances, "-o", label='Structural change')
    plt.title('Distance from Initial Structure')
    plt.xlabel('Optimization Step')
    plt.ylabel('RMSD from Initial Structure')
    plt.xlim(0, len(distances) + 1)
    plt.legend()
    plt.savefig('conv_structure.png', bbox_inches='tight')

    print("处理完成！")

if __name__ == "__main__":
    # 检查输入参数
    if len(sys.argv) < 2:
        print("用法: python convergence_plot.py ./opt/ (POSCAR_ini)")
        sys.exit(1)

    # 读取命令行参数
    path = sys.argv[1]
    initial_structure_path = sys.argv[2] if len(sys.argv) == 3 else None

    # 绘制收敛图
    plot_convergence(path, initial_structure_path)
