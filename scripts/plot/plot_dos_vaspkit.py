"""
文件名: plot_dos_vaspkit.py
运行方法: python plot_dos_vaspkit.py Si Sb
功能描述: 该脚本读取 TDOS 和 PDOS 数据文件，并绘制总态密度 (TDOS) 和部分态密度 (PDOS) 图。
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
import sys

# 设置 matplotlib 的全局绘图参数
plt.rcParams['font.size'] = 14
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'

# 获取命令行参数中指定的元素列表
elements = sys.argv[1:]

# 读取 TDOS 数据文件
with open('./TDOS.dat', 'r') as f:
    data = []
    for line in f.readlines():
        l = []
        # 去除换行符并分割每一行的内容
        for item in line.replace("\n", "").split(" "):
            if item != "":
                l.append(item)
        data.append(l)
    # 将数据转换为 Pandas DataFrame
    df0 = pd.DataFrame(data)

# 将数据类型转换为数值型
x0 = pd.to_numeric(df0[0][1:].reset_index(drop=True), errors='coerce')
y0 = pd.to_numeric(df0[1][1:].reset_index(drop=True), errors='coerce')

# 绘制 TDOS 图像
plt.plot(x0, y0, c="k", linewidth=1, label="TDOS")
plt.fill_between(x0, 0, y0, facecolor='grey', alpha=0.3)

def plot_pdos(path):
    """
    读取 PDOS 数据并绘图。

    参数:
    path (str): PDOS 文件的路径
    """
    # 获取文件名（不含扩展名）作为标签
    label = os.path.splitext(path)[0]
    # 读取 PDOS 数据文件
    with open(path, 'r') as f:
        data = []
        for line in f.readlines():
            l = []
            for item in line.replace("\n", "").split(" "):
                if item != "":
                    l.append(item)
            data.append(l)
        df1 = pd.DataFrame(data)
    
    # 提取能量和 DOS 数据
    x1 = pd.to_numeric(df1[0][1:].reset_index(drop=True), errors='coerce')
    y11 = np.zeros_like(x1)
    
    # 将 4 列 DOS 数据相加
    for i in range(4):
        y1 = pd.to_numeric(df1[i+1][1:].reset_index(drop=True), errors='coerce')
        y11 += y1

    # 绘制 PDOS 图像
    plt.plot(x1, y11, alpha=0.7, linewidth=1, label=label)

# 为每个指定的元素绘制 PDOS 图像
for element in elements:
    path = f'PDOS_{element}.dat'
    plot_pdos(path)

# 设置图形标签和范围
plt.xlabel('Energy (eV)')
plt.ylabel('DOS')
plt.xlim(x0.min(), x0.max())
plt.ylim(0, y0.max() * 1.1)

# 显示图例
plt.legend(loc='best')

# 保存图像为 PNG 文件
plt.savefig('dos.png', bbox_inches='tight')
