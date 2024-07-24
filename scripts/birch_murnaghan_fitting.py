#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
用法: 运行脚本前，请确保当前目录下有ev.csv文件。

该脚本从ev.csv文件中读取晶格常数和能量数据，并使用Birch-Murnaghan方程进行拟合。拟合结果将生成一个散点图和拟合曲线，并保存为ev.png。
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import leastsq

# 设置画图格式
plt.rcParams['font.size'] = 14
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'

def Birch_Murnaghan(p, x):
    """
    Birch-Murnaghan状态方程

    参数:
        p: 参数列表 [E0, V0, B0, B1]
        x: 晶格常数数组

    返回:
        对应的能量值数组
    """
    E0, V0, B0, B1 = p
    eta = (V0 / (x ** 3)) ** (2 / 3)
    return E0 + (9 * V0 * B0 / 16) * (((eta - 1) ** 3) * B1 + ((eta - 1) ** 2) * (6 - 4 * eta))

def error(p, x, y):
    """
    误差函数

    参数:
        p: 参数列表 [E0, V0, B0, B1]
        x: 晶格常数数组
        y: 实验能量值数组

    返回:
        拟合误差
    """
    return Birch_Murnaghan(p, x) - y

# 读取数据
f = "./ev.csv"
df = pd.read_csv(f)
x = np.array(df['lat'])
y = np.array(df['ene'])

# 初始参数猜测
p0 = [min(y), min(x) ** 3, 100, 1]

# 最小二乘法拟合
para = leastsq(error, p0, args=(x, y))
y_fitted = Birch_Murnaghan(para[0], x)

# 输出拟合结果
print("Emin =", para[0][0])
print("opt_lat =", para[0][1] ** (1 / 3))
print("B0 =", para[0][2])
print("B1 =", para[0][3])

# 绘图
plt.figure(figsize=(8, 6), dpi=300)
plt.scatter(x, y, 20, c='r', label="DFT")
plt.plot(x, y_fitted, '-b', label="B-M Fitting")
plt.xlabel("Lattice constant", fontsize=20)
plt.ylabel("Energy (eV)", fontsize=20)
plt.legend()
plt.savefig('./ev.png', bbox_inches='tight')
