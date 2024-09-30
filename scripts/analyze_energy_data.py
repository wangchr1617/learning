"""
文件名: analyze_energy_data.py
运行方法: python analyze_energy_data.py
功能描述: 读取 NEP 和 DFT 能量数据，计算每个结构的标准差，按标准差降序排列，选择和绘制指定范围内的数据。
注意：结合脚本 nep_evaluate.py 使用。
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from ase.io import read, write
import warnings

# 设置 matplotlib 全局参数
plt.rcParams['font.size'] = 14
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'

# 忽略警告
warnings.filterwarnings('ignore')

def load_data():
    """
    读取 NEP 和 DFT 能量数据，计算标准差，并按标准差降序排列。
    """
    nep = np.loadtxt('./nep_energy-1.txt', usecols=(1)).reshape(-1, 1)
    dft = np.loadtxt('./dft_energy.txt', usecols=(0)).reshape(-1, 1)
    natoms = np.loadtxt('./atoms_num.txt', usecols=(0)).reshape(-1, 1)
    data = np.concatenate((nep, dft, natoms), axis=1)
    
    df = pd.DataFrame(data, columns=["nep", "dft", "num"])
    with open('./config_type.txt', 'r') as file:
        df["config_type"] = file.read().splitlines()[:len(df)]
    
    df["dft"] = df["dft"] / df["num"]
    df["std"] = df[["nep", "dft"]].std(axis=1)
    df.sort_values(by="std", ascending=False, inplace=True)
    df.to_csv('./std.csv', index=False)
    
    return df

def select_data(df, begin, end=None):
    """
    选择指定范围的数据，并保存到文件。
    """
    selected_indices = df.index.tolist()[begin:end]
    structures = read('./model.xyz', ':')
    write('selected.xyz', [structures[i] for i in selected_indices])

def plot_energy(df, begin=None, end=None):
    """
    绘制 DFT 和 NEP 能量比较图。
    """
    plt.figure(figsize=(5, 4), dpi=300)
    x = df['dft'][begin:end]
    y = df['nep'][begin:end]
    xmin, xmax = min(x.min(), y.min()) * 1.005, max(x.max(), y.max()) * 0.995
    
    plt.plot([xmin, xmax], [xmin, xmax], linewidth=1, color='gray')
    plt.scatter(x, y, marker='o', s=1.5)
    plt.xlabel('DFT energy (eV/atom)')
    plt.ylabel('NEP energy (eV/atom)')
    plt.xlim(xmin, xmax)
    plt.ylim(xmin, xmax)
    plt.tight_layout()
    plt.savefig('./energy.png', bbox_inches='tight')

def plot_energy_population(df):
    """
    绘制能量密度分布图。
    """
    dft = np.round(df['dft'], decimals=3)
    unique_energies, counts = np.unique(dft, return_counts=True)

    print(f"DFT Energy: max = {unique_energies.max()}, min = {unique_energies.min()}")
    
    plt.figure(figsize=(5, 3))
    plt.plot(unique_energies, counts, color='r', label='DFT')
    plt.xlabel('Energy per atom (eV/atom)')
    plt.ylabel('Density of population')
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig('./energy_population.png', bbox_inches='tight')

def main():
    df = load_data()
    plot_energy(df)
    plot_energy_population(df)
    # select_data(df, 5)

if __name__ == '__main__':
    main()
