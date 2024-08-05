"""
用法: python rag_energy.py

功能描述:
该脚本从 `dft_energy.txt` 和 `atoms_num.txt` 文件中读取数据，计算并绘制能量密度分布图。
图表中将比较未松弛和松弛后的能量分布。

输入文件:
- `dft_energy.txt`: 包含每个模型的能量数据。
- `atoms_num.txt`: 包含每个模型的原子数量。

结果:
- 生成能量密度分布图，并保存为 `energy_population.png`。
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# 设置图形参数
plt.rcParams['font.size'] = 14
plt.rcParams["axes.unicode_minus"] = False
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'

import warnings
warnings.filterwarnings('ignore')

def qbc():
    """
    读取能量和原子数量数据，计算每个原子的能量。

    返回:
    - df: 包含能量和原子数量的数据框，筛选出特定原子数量的行。
    """
    # 读取数据
    ene = np.loadtxt('./dft_energy.txt', usecols=(0)).reshape((-1, 1))
    natoms = np.loadtxt('./atoms_num.txt', usecols=(0)).reshape((-1, 1))
    
    # 将能量和原子数量合并为一个数组
    data = np.concatenate((ene, natoms), axis=1)
    
    # 创建DataFrame并筛选原子数量为62的行
    df = pd.DataFrame(data, columns=["ene", "natoms"])
    df = df[df["natoms"] == 62]
    
    # 计算每个原子的能量
    df["ene"] = np.array(df["ene"]) / np.array(df["natoms"])
    
    return df

def count(a, x):
    """
    计算数组a中每个x值出现的次数。

    参数:
    - a: 要计数的数组。
    - x: 要计数的值的列表。

    返回:
    - 计数结果的数组。
    """
    return np.array([np.sum(a == i) for i in x])

def plot(df, step=20):
    """
    绘制能量密度分布图。

    参数:
    - df: 包含能量和原子数量的数据框。
    - step: 用于计算能量密度的步长。
    """
    # 计算初始能量和松弛能量
    ini = np.round(np.array(df['ene'])[::step], decimals=3)
    rag = np.round(np.array(df['ene'])[::5], decimals=3)
    
    # 获取不重复的能量值并计数
    x1 = np.unique(sorted(ini))
    y1 = count(ini, x1)
    
    x2 = np.unique(sorted(rag))
    y2 = count(rag, x2)
    
    # 输出最大最小能量
    print(f"Initial Energy: max = {max(x1)}, min = {min(x1)}")
    print(f"Relaxed Energy: max = {max(x2)}, min = {min(x2)}")
    
    # 绘图
    plt.figure(figsize=(5, 3))
    plt.plot(x1, y1, color='black', label='Not relaxed')
    plt.plot(x2, y2, color='red', label='Relaxed')
    plt.xlim(-4.25, -4.18)
    plt.xlabel('Energy per atom (eV/atom)')
    plt.ylabel('Density of population')
    plt.legend()
    plt.grid(True)
    plt.savefig('./energy_population.png', bbox_inches='tight')
    print("图像已保存为 'energy_population.png'")

def main():
    df = qbc()
    plot(df)

if __name__ == "__main__":
    main()
