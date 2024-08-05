# 文件名: analyze_energy_data.py
# 运行方法: python analyze_energy_data.py
# 功能描述: 读取NEP和DFT能量数据，计算每个结构的标准差，按标准差降序排列，并选择和绘制指定范围内的数据。

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from ase.io import read, write
import warnings

# 设置matplotlib全局参数
plt.rcParams['font.size'] = 14
plt.rcParams["axes.unicode_minus"] = False
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'

# 忽略警告
warnings.filterwarnings('ignore')

def qbc():
    """
    读取NEP和DFT能量数据，计算每个结构的标准差，并按标准差降序排列。
    
    返回:
        pd.DataFrame: 包含能量数据和标准差的数据框。
    """
    n1 = np.loadtxt('./nep_energy-1.txt', usecols=(1)).reshape((-1, 1))
    n4 = np.loadtxt('./dft_energy.txt', usecols=(0)).reshape((-1, 1))
    n5 = np.loadtxt('./atoms_num.txt', usecols=(0)).reshape((-1, 1))
    n = np.concatenate((n1, n4, n5), axis=1)
    df = pd.DataFrame(n, columns=["nep1", "dft", "num"])
    with open('./config_type.txt', 'r') as file:
        config_type = file.read().splitlines()
    df["config_type"] = config_type[:len(df)]
    df["dft"] = np.array(df["dft"]) / np.array(df["num"])
    df["std"] = df[["nep1", "dft"]].std(axis=1, numeric_only=True)
    df = df.sort_values(by=["std"], ascending=[False])
    df.to_csv('./std.csv', index=False)
    return df

def read_df():
    """
    读取保存的标准差数据文件。
    
    返回:
        pd.DataFrame: 标准差数据框。
    """
    df = pd.read_csv('std.csv')
    return df

def select(df, begin, end=None):
    """
    从标准差数据框中选择指定范围内的数据，并保存到新文件。
    
    参数:
        df (pd.DataFrame): 标准差数据框。
        begin (int): 开始索引。
        end (int, 可选): 结束索引，默认为None。
    """
    df_index_list = df.index.tolist()[begin:end]
    a = read('./model.xyz', ':')
    write('qbc.xyz', [a[i] for i in df_index_list])

def plot(df, begin=None, end=None):
    """
    绘制DFT和NEP能量比较图。
    
    参数:
        df (pd.DataFrame): 标准差数据框。
        begin (int, 可选): 开始索引，默认为None。
        end (int, 可选): 结束索引，默认为None。
    """
    plt.figure(figsize=(5, 4), dpi=300)
    x = df['dft'][begin:end]
    y = df['nep1'][begin:end]
    xmax = max([max(x), max(y)]) * 0.995
    xmin = min([min(x), min(y)]) * 1.005
    plt.plot([xmin, xmax], [xmin, xmax], linewidth=1)
    plt.scatter(x, y, marker='o', s=1.5)
    plt.xlabel('DFT energy (eV/atom)')
    plt.ylabel('NEP energy (eV/atom)')
    plt.xlim(xmin, xmax)
    plt.ylim(xmin, xmax)
    plt.tight_layout()
    plt.savefig('./energy.png', bbox_inches='tight')

def main():
    """
    主函数，调用其他函数实现数据处理和绘图。
    """
    df = qbc()
    # 如果需要读取已有的标准差数据文件，可以注释上行，取消下行注释
    # df = read_df()
    # 如果需要选择部分数据，可以取消下行注释并设置范围
    # select(df, begin=50)
    plot(df)

if __name__ == '__main__':
    main()
