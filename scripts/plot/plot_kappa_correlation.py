"""
文件名: plot_kappa_correlation.py
运行方法: 直接在支持 Python 环境下运行此脚本
功能描述: 该脚本加载 HAC 数据并计算平均值，然后绘制 $\kappa^{in}$、$\kappa^{out}$ 和总 $\kappa$ 的相关性图。
"""

import numpy as np
import matplotlib.pyplot as plt
from gpyumd.load import load_hac

# 设置画图参数
plt.rcParams['font.size'] = 14
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'

# 加载数据
hac = load_hac([50000] * 3, [10] * 3)

# 打印运行信息
print("Runs:", list(hac.keys()))
print("Run Data:", list(hac['run0'].keys()))

# 提取时间数据
t = hac['run0']['t']

# 初始化数据平均值数组
hac_ave_i = np.zeros(hac['run0']['jxijx'].shape[0])
hac_ave_o = np.zeros_like(hac_ave_i)
ki_ave, ko_ave = np.zeros_like(hac_ave_i), np.zeros_like(hac_ave_o)

# 计算每次运行的平均值
for runkey in hac.keys():
    hac_ave_i += hac[runkey]['jxijx'] + hac[runkey]['jyijy']
    hac_ave_o += hac[runkey]['jxojx'] + hac[runkey]['jyojy']
    ki_ave += hac[runkey]['kxi'] + hac[runkey]['kyi']
    ko_ave += hac[runkey]['kxo'] + hac[runkey]['kyo']

# 标准化处理
hac_ave_i /= hac_ave_i.max()
hac_ave_o /= hac_ave_o.max()
ki_ave /= 6.0
ko_ave /= 6.0

# 设置字体和线宽
aw = 2
fs = 16
font = {'size': fs}
plt.rc('font', **font)
plt.rc('axes', linewidth=aw)

# 设置图形属性
def set_fig_properties(ax_list):
    tl = 8
    tw = 2
    tlm = 4

    for ax in ax_list:
        ax.tick_params(which='major', length=tl, width=tw)
        ax.tick_params(which='minor', length=tlm, width=tw)
        ax.tick_params(which='both', axis='both', direction='in', right=True, top=True)

# 创建图形并绘图
plt.figure(figsize=(12, 4))

# 子图1：kappa^in
plt.subplot(1, 3, 1)
set_fig_properties([plt.gca()])
for runkey in hac.keys():
    plt.plot(hac[runkey]['t'], (hac[runkey]['kxi'] + hac[runkey]['kyi']) / 2, color='C7', alpha=0.5)
plt.plot(t, ki_ave, color='C3', linewidth=1)
plt.xlim([0, 1000])
plt.gca().set_xticks(range(0, 1001, 200))
plt.ylim([0, 5])
plt.xlabel('Correlation Time (ps)')
plt.ylabel(r'$\kappa^{in}$ (W/m/K)')
plt.title('(a)')

# 子图2：kappa^out
plt.subplot(1, 3, 2)
set_fig_properties([plt.gca()])
for runkey in hac.keys():
    plt.plot(hac[runkey]['t'], (hac[runkey]['kxo'] + hac[runkey]['kyo']) / 2, color='C7', alpha=0.5)
plt.plot(t, ko_ave, color='C0', linewidth=1)
plt.xlim([0, 1000])
plt.gca().set_xticks(range(0, 1001, 200))
plt.ylim([0, 5])
plt.xlabel('Correlation Time (ps)')
plt.ylabel(r'$\kappa^{out}$ (W/m/K)')
plt.title('(b)')

# 子图3：kappa 总和
plt.subplot(1, 3, 3)
set_fig_properties([plt.gca()])
plt.plot(t, ko_ave, color='C0', linewidth=1)
plt.plot(t, ki_ave, color='C3', linewidth=1)
plt.plot(t, ki_ave + ko_ave, color='k', linewidth=1)
plt.xlim([0, 1000])
plt.gca().set_xticks(range(0, 1001, 200))
plt.ylim([0, 5])
plt.xlabel('Correlation Time (ps)')
plt.ylabel(r'$\kappa$ (W/m/K)')
plt.title('(c)')

# 调整布局并保存图像
plt.tight_layout()
plt.savefig('./EMD.png', bbox_inches='tight')
