"""
文件名：plot_dos_gpumd.py
运行方法：python plot_dos_gpumd.py
功能描述：加载并可视化 GPUMD 输出的声子态密度 (DOS) 和速度自相关函数 (VAC) 的数据，生成对应的图像并保存。
"""

import numpy as np
import matplotlib.pyplot as plt
from gpyumd.load import load_dos, load_vac

# 设置全局图形参数
plt.rcParams['font.size'] = 14
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'

def set_fig_properties(ax_list):
    """
    设置子图的刻度属性。

    参数：
    ax_list (list): Matplotlib的Axes对象列表。
    """
    tick_length_major = 8  # 主刻度线长度
    tick_width = 2         # 刻度线宽度
    tick_length_minor = 4  # 次刻度线长度

    for ax in ax_list:
        ax.tick_params(which='major', length=tick_length_major, width=tick_width)
        ax.tick_params(which='minor', length=tick_length_minor, width=tick_width)
        ax.tick_params(which='both', axis='both', direction='in', right=True, top=True)

# 加载DOS和VAC数据
num_corr_steps = 1000
dos = load_dos(num_dos_points=num_corr_steps)['run0']
vac = load_vac(num_corr_steps)['run0']

# 计算总的DOS和VAC
dos['DOSxyz'] = dos['DOSx'] + dos['DOSy'] + dos['DOSz']
vac['VACxyz'] = vac['VACx'] + vac['VACy'] + vac['VACz']
vac['VACxyz'] /= vac['VACxyz'].max()  # 归一化处理

# 输出DOS和VAC的键值
print('DOS:', list(dos.keys()))
print('VAC:', list(vac.keys()))

# 设置画布大小
plt.figure(figsize=(12, 10))

# 绘制子图1：VAC的x, y, z方向分量
plt.subplot(2, 2, 1)
plt.plot(vac['t'], vac['VACx'] / vac['VACx'].max(), color='C3', linewidth=3)
plt.plot(vac['t'], vac['VACy'] / vac['VACy'].max(), color='C0', linestyle='--', linewidth=3)
plt.plot(vac['t'], vac['VACz'] / vac['VACz'].max(), color='C2', linestyle='-.', zorder=100, linewidth=3)
plt.xlabel('Correlation Time (ps)')
plt.ylabel('VAC (Normalized)')
plt.xlim(0, 5)
plt.legend(['x', 'y', 'z'])
plt.title('(a)')

# 绘制子图2：DOS的x, y, z方向分量
plt.subplot(2, 2, 2)
plt.plot(dos['nu'], dos['DOSx'], color='C3', linewidth=3)
plt.plot(dos['nu'], dos['DOSy'], color='C0', linestyle='--', linewidth=3)
plt.plot(dos['nu'], dos['DOSz'], color='C2', linestyle='-.', zorder=100, linewidth=3)
plt.xlabel(r'$\nu$ (THz)')
plt.ylabel('PDOS (1/THz)')
plt.xlim(0, 8)
plt.ylim(0, None)
plt.legend(['x', 'y', 'z'])
plt.title('(b)')

# 绘制子图3：总的VAC
plt.subplot(2, 2, 3)
plt.plot(vac['t'], vac['VACxyz'], color='k', linewidth=3)
plt.xlabel('Correlation Time (ps)')
plt.ylabel('VAC (Normalized)')
plt.xlim(0, 5)
plt.title('(c)')

# 绘制子图4：总的DOS
plt.subplot(2, 2, 4)
plt.plot(dos['nu'], dos['DOSxyz'], color='k', linewidth=3)
plt.xlabel(r'$\nu$ (THz)')
plt.ylabel('PDOS (1/THz)')
plt.xlim(0, 8)
plt.ylim(0, None)
plt.title('(d)')

# 保存图像
plt.tight_layout()
plt.savefig('./VDOS.png', bbox_inches='tight')
