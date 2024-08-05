# 文件名: nep_evaluate.py
# 运行方法: python nep_evaluate.py nep_1.txt model.xyz 1
# 功能描述: 使用NEP计算器计算给定模型的能量，并将结果保存到文件中。

from ase.io import read
from pynep.calculate import NEP
import numpy as np
import sys

def evaluate(c, t, label):
    """
    使用NEP计算器计算模型的能量。

    参数:
        c (str): NEP配置文件的路径。
        t (str): 模型文件的路径。
        label (str): 输出文件的标签。
    """
    # 初始化NEP计算器
    calc = NEP(c)
    print(calc)

    # 读取模型轨迹
    traj = read(t, ':')

    # 初始化能量和力的列表
    e, f = [], []
    idx = 0

    # 遍历每个原子结构，计算能量和力
    for atoms in traj:
        num = len(atoms)  # 获取原子数目
        atoms.set_calculator(calc)  # 设置计算器
        e.append([idx, atoms.get_potential_energy() / num])  # 计算每个原子的能量
        # f.append(atoms.get_forces().reshape(-1))  # 计算力并拉平为一维数组
        idx += 1

    # 转换能量和力为numpy数组
    e = np.array(e)
    # f = np.concatenate(f)

    # 保存计算结果到文件
    np.savetxt('nep_energy-{}.txt'.format(label), e, delimiter=' ', fmt="%.4f")
    # np.savetxt('nep_force-{}.txt'.format(label), f, delimiter=' ', fmt="%.4f")

def main():
    """
    主函数，解析命令行参数并调用评估函数。
    """
    if len(sys.argv) != 4:
        print("Usage: python nep_evaluate.py <nep_file> <model_file> <label>")
        sys.exit(1)

    nep_file = sys.argv[1]
    model_file = sys.argv[2]
    label = sys.argv[3]

    evaluate(nep_file, model_file, label)

if __name__ == '__main__':
    main()
