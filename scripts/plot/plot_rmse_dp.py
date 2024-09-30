"""
文件名: plot_rmse_dp.py
运行方法: 
dp test -m frozen_model.pb -s ./train/Ge64Te64/ -n 100 -d results
python plot_rmse_dp.py
功能描述: 
该脚本用于绘制 DFT 与 DP 计算结果的能量和力的比较图。
"""

import numpy as np
import matplotlib.pyplot as plt

# 统一画图格式
plt.rcParams['font.size'] = 14
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'

def plot(ax, data, key, xlabel, ylabel, min_val, max_val):
    """
    绘制散点图，并在对角线上绘制 y=x 的红线。
    
    参数:
        ax: Matplotlib 的轴对象
        data: 包含数据的字典
        key: 字符串，指示绘图数据的键
        xlabel: x 轴标签
        ylabel: y 轴标签
        min_val: x 和 y 的最小值，用于设置坐标轴范围
        max_val: x 和 y 的最大值，用于设置坐标轴范围
    """
    data_key = f'data_{key}'
    pred_key = f'pred_{key}'
    ax.scatter(data[data_key], data[pred_key], label=key, s=6)
    ax.legend()
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_xlim(min_val, max_val)
    ax.set_ylim(min_val, max_val)
    ax.plot([min_val, max_val], [min_val, max_val], 'r', lw=1)

# 原子数，注意根据实际情况修改
natom = 128

# 读取数据文件
data_e = np.genfromtxt("results.e.out", names=["data_e", "pred_e"])
data_f = np.genfromtxt("results.f.out", names=["data_fx", "data_fy", "data_fz", "pred_fx", "pred_fy", "pred_fz"])

# 单位换算：将能量除以原子数
data_e['data_e'] /= natom
data_e['pred_e'] /= natom

# 计算能量和力的最小值和最大值，用于设定坐标轴范围
data_e_stacked = np.column_stack((data_e['data_e'], data_e['pred_e']))
data_f_stacked = np.column_stack((data_f['data_fx'], data_f['data_fy'], data_f['data_fz'], data_f['pred_fx'], data_f['pred_fy'], data_f['pred_fz']))

min_val_e, max_val_e = np.min(data_e_stacked), np.max(data_e_stacked)
min_val_f, max_val_f = np.min(data_f_stacked), np.max(data_f_stacked)

# 创建画布并绘制图像
fig, axs = plt.subplots(1, 2, figsize=(12, 5))
plot(axs[0], data_e, 'e', 'DFT energy (eV/atom)', 'DP energy (eV/atom)', min_val_e, max_val_e)

for force_direction in ['fx', 'fy', 'fz']:
    plot(axs[1], data_f, force_direction, 'DFT force (eV/Å)', 'DP force (eV/Å)', min_val_f, max_val_f)

# 保存图片为 PNG 文件
plt.savefig('DP_vs_DFT.png', bbox_inches='tight')
