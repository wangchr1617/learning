# 文件名: plot_force_error.py
# 运行方法: python plot_force_error.py
# 功能描述: 计算参考数据、DFT数据和NEP数据之间的力误差，并绘制误差的直方图和拟合的正态分布曲线。

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import norm

# 设置matplotlib全局参数
plt.rcParams['font.size'] = 14
plt.rcParams["axes.unicode_minus"] = False
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'

import warnings
warnings.filterwarnings('ignore')

def calculate_force_error():
    """
    计算参考数据、DFT数据和NEP数据之间的力误差。

    返回:
        DataFrame: 包含力误差数据的DataFrame。
    """
    # 加载参考、DFT和NEP数据
    r = np.loadtxt('./ref.txt').reshape((-1, 1))
    d = np.loadtxt('./dft.txt').reshape((-1, 1))
    n = np.loadtxt('./nep.txt').reshape((-1, 1))

    # 将数据合并为DataFrame
    n = np.concatenate((r, d, n), axis=1)
    df = pd.DataFrame(n, columns=["ref", "dft", "nep"])

    # 计算力误差
    df["e_d"] = df["ref"] - df["dft"]
    df["e_n"] = df["ref"] - df["nep"]

    return df

def plot_force_error(df):
    """
    绘制力误差的直方图和拟合的正态分布曲线。

    参数:
        df (DataFrame): 包含力误差数据的DataFrame。
    """
    # 提取力误差数据
    e_d = df['e_d'].to_numpy()
    e_n = df['e_n'].to_numpy()
   
    plt.figure(figsize=(5, 4), dpi=300)

    # 绘制DFT力误差的直方图和拟合曲线
    plt.hist(e_d, bins=30, color='k', alpha=0.3, density=True, label='DFT')
    xmin, xmax = plt.xlim()
    x_d = np.linspace(xmin, xmax, 100)
    mu_d, std_d = norm.fit(e_d)
    p_d = norm.pdf(x_d, mu_d, std_d)
    plt.plot(x_d, p_d, linewidth=1, linestyle='--', color='k', alpha=0.75)
    plt.fill_between(x_d, p_d, where=(p_d >= 0), color='k', alpha=0.5)    

    # 绘制NEP力误差的直方图和拟合曲线
    plt.hist(e_n, bins=30, color='r', alpha=0.3, density=True, label='NEP')
    xmin, xmax = plt.xlim()
    x_n = np.linspace(xmin, xmax, 100)
    mu_n, std_n = norm.fit(e_n)
    p_n = norm.pdf(x_n, mu_n, std_n)
    plt.plot(x_n, p_n, linewidth=1, linestyle='--', color='r', alpha=0.75)
    plt.fill_between(x_n, p_n, where=(p_n >= 0), color='r', alpha=0.5)    

    # 设置图表的标题和标签
    plt.xlabel('Force error (eV/A)')
    plt.ylabel('Density')
    plt.xlim(-0.25, 0.25)
    plt.xticks(np.arange(-0.2, 0.25, 0.1))
    plt.ylim(0, max(p_d))
    plt.legend()
    
    plt.tight_layout()
    plt.savefig('./Force_error.png', bbox_inches='tight')
    print("Force error plot has been saved as Force_error.png.")

def main():
    """
    主函数，计算力误差并绘图。
    """
    df = calculate_force_error()
    plot_force_error(df)

if __name__ == '__main__':
    main()
