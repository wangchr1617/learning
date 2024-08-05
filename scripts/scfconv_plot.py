"""
用法: python scfconv_plot.py [step_index]

功能描述:
该脚本用于从VASP的OSZICAR文件中读取自洽步骤的能量变化，并绘制能量收敛图。

参数说明:
- [step_index]: 要绘制的自洽步骤索引。

结果:
- 生成能量收敛图，并保存为 `scf.png` 文件。
"""

from ase.io import read
import numpy as np
import matplotlib.pyplot as plt
import sys
import warnings

# 忽略警告
warnings.filterwarnings('ignore')

# 设置图形参数
plt.rcParams['font.size'] = 14
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'

def main():
    # 检查命令行参数
    if len(sys.argv) != 2:
        print("用法: python scfconv_plot.py [step_index]")
        sys.exit(1)

    # 获取自洽步骤索引
    idx = int(sys.argv[1])

    # 初始化自洽步骤和能量列表
    sc_steps_list = []
    sc_energies_list = []

    # 读取OSZICAR文件，提取自洽步骤和能量
    with open("OSZICAR", "r") as f:
        i = 0
        for line in f.readlines():
            if i > idx:
                break
            if "N " in line:
                # 初始化新的步骤和能量列表
                sc_steps = []
                sc_energies = []
            elif "E0=" in line:
                # 存储当前步骤的能量数据
                sc_steps_list.append(sc_steps)
                sc_energies_list.append(sc_energies)
                i += 1
            else:
                # 提取步骤和能量数据
                parts = line.split()
                sc_steps.append(int(parts[1]))
                sc_energies.append(float(parts[2]))

    # 绘制能量收敛图
    plt.figure(figsize=(10, 6))
    plt.grid(True)
    plt.plot(sc_steps_list[idx], sc_energies_list[idx], '-o', label=f'Step {idx}')
    plt.title('Energy of each self-consistency step')
    plt.xlabel('Self-consistency steps')
    plt.ylabel('Energy (eV)')
    plt.xlim(0, max(sc_steps_list[idx]) + 1)
    plt.legend()
    plt.savefig('scf.png', bbox_inches='tight')
    print("图像已保存为 'scf.png'.")

if __name__ == "__main__":
    main()
