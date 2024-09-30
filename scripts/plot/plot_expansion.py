'''
文件名: plot_lattice_expansion.py
运行方法: python plot_lattice_expansion.py
功能描述: 根据输入的热力学数据文件生成晶格膨胀的相关图表，并保存为 PNG 图片。
'''

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# 设置画图参数
plt.rcParams['font.size'] = 14
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'

# 配色方案
palette = [
    [0.9678, 0.4413, 0.5358],
    [0.8088, 0.5635, 0.1950],
    [0.5921, 0.6418, 0.1935],
    [0.1978, 0.6956, 0.3995],
    [0.2104, 0.6773, 0.6434],
    [0.2234, 0.6566, 0.8171],
    [0.6423, 0.5498, 0.9583],
    [0.9604, 0.3814, 0.8683]
]

# 计算体积
def volume(a, b, c, ncell):
    arr = np.array([a, b, c]).T
    return np.linalg.det(arr) / ncell

# 计算晶格长度
def length(vec, ncell):
    return np.linalg.norm(vec) / ncell

# 计算晶格角度
def angle(vec1, vec2):
    cos_theta = np.dot(vec1, vec2) / (np.linalg.norm(vec1) * np.linalg.norm(vec2))
    angle_radians = np.arccos(np.clip(cos_theta, -1.0, 1.0))
    return np.abs(90 - np.degrees(angle_radians))

# 生成并保存图表
def plot_md(filename):
    # 读取数据文件
    data = np.loadtxt(filename)
    cols = ["T", "K", "U", "Px", "Py", "Pz", "Pyz", "Pxz", "Pxy", 
            "ax", "ay", "az", "bx", "by", "bz", "cx", "cy", "cz"]
    thermo = pd.DataFrame(data, columns=cols)
    
    ncell = 12
    natom = (ncell**3) * 8
    temp = np.arange(100, 800, 1)[::-1]
    interval = int(len(thermo) / len(temp))
    thermo = thermo[::interval].reset_index(drop=True)
    
    # 计算体积和晶格参数
    v_ave = np.array([
        volume(thermo.loc[i, ['ax', 'ay', 'az']], 
               thermo.loc[i, ['bx', 'by', 'bz']], 
               thermo.loc[i, ['cx', 'cy', 'cz']], ncell**3)
        for i in range(len(thermo))
    ])
    
    a_ave = np.array([length(thermo.loc[i, ['ax', 'ay', 'az']], ncell) for i in range(len(thermo))])
    b_ave = np.array([length(thermo.loc[i, ['bx', 'by', 'bz']], ncell) for i in range(len(thermo))])
    c_ave = np.array([length(thermo.loc[i, ['cx', 'cy', 'cz']], ncell) for i in range(len(thermo))])
    
    alpha = np.array([angle(thermo.loc[i, ['bx', 'by', 'bz']], thermo.loc[i, ['cx', 'cy', 'cz']]) 
                      for i in range(len(thermo))])
    beta = np.array([angle(thermo.loc[i, ['ax', 'ay', 'az']], thermo.loc[i, ['cx', 'cy', 'cz']]) 
                     for i in range(len(thermo))])
    gamma = np.array([angle(thermo.loc[i, ['ax', 'ay', 'az']], thermo.loc[i, ['bx', 'by', 'bz']]) 
                      for i in range(len(thermo))])
    
    p_ave = thermo[['Px', 'Py', 'Pz', 'Pyz', 'Pxz', 'Pxy']].mean(axis=1)
    tot_ene = (thermo['K'] + thermo['U']) / natom
    
    # 创建图表
    plt.figure(figsize=(10, 8), dpi=300)
    
    plt.subplot(2, 2, 1)
    plt.plot(temp, v_ave, c='k')
    plt.axvline(275, c='k', linestyle='--', linewidth=1, alpha=0.5)
    plt.axvline(450, c='k', linestyle='--', linewidth=1, alpha=0.5)
    plt.xlim(100, 800)
    plt.xlabel('Temperature (K)')
    plt.ylabel(r'Normalized Lattice Volume (Å³/f.u.)')
    plt.xticks(np.arange(100, 850, 100))
    plt.title('(a)')
    
    plt.subplot(2, 2, 2)
    plt.plot(temp, a_ave, c=palette[1], label='a')
    plt.plot(temp, b_ave, c=palette[2], label='b')
    plt.plot(temp, c_ave, c=palette[3], label='c')
    plt.axvline(275, c='k', linestyle='--', linewidth=1, alpha=0.5)
    plt.axvline(450, c='k', linestyle='--', linewidth=1, alpha=0.5)
    plt.xlim(100, 800)
    plt.xlabel('Temperature (K)')
    plt.ylabel(r'Normalized Lattice Length (Å)')
    plt.xticks(np.arange(100, 850, 100))
    plt.legend(loc="lower right")
    plt.title('(b)')
    
    plt.subplot(2, 2, 3)
    plt.plot(temp, alpha, c=palette[1], label=r'$\alpha$')
    plt.plot(temp, beta,  c=palette[2], label=r'$\beta$')
    plt.plot(temp, gamma, c=palette[3], label=r'$\gamma$')
    plt.axvline(275, c='k', linestyle='--', linewidth=1, alpha=0.5)
    plt.axvline(450, c='k', linestyle='--', linewidth=1, alpha=0.5)
    plt.xlim(100, 800)
    plt.ylim(-0.1, )
    plt.xlabel('Temperature (K)')
    plt.ylabel(r'90°- Lattice Angle (°)')
    plt.xticks(np.arange(100, 850, 100))
    plt.legend(loc="upper right")
    plt.title('(c)')
    
    plt.subplot(2, 2, 4)
    plt.plot(temp, thermo['Px'], c=palette[0], alpha=0.5, label='Px')
    plt.plot(temp, thermo['Py'], c=palette[1], alpha=0.5, label='Py')
    plt.plot(temp, thermo['Pz'], c=palette[2], alpha=0.5, label='Pz')
    plt.plot(temp, thermo['Pyz'], c=palette[3], alpha=0.5, label='Pyz')
    plt.plot(temp, thermo['Pxz'], c=palette[4], alpha=0.5, label='Pxz')
    plt.plot(temp, thermo['Pxy'], c=palette[5], alpha=0.5, label='Pxy')
    plt.plot(temp, p_ave, c='k', label='P_ave')
    plt.xlim(100, 800)
    plt.ylim(-0.5, 0.5)
    plt.xlabel('Temperature (K)')
    plt.ylabel('Pressure (GPa)')
    plt.xticks(np.arange(100, 850, 100))
    plt.legend(loc="upper right")
    plt.title('(d)')
    
    plt.tight_layout()
    # 保存图表为 PNG 格式
    plt.savefig('./lattice_expansion.png', bbox_inches='tight')

# 调用函数，生成图表
plot_md('./thermo.out')
