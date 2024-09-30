"""
文件名: plot_loss_dp.py
运行方法: python plot_loss_dp.py
功能描述: 该脚本用于绘制 DeePMD 训练过程中的损失曲线。
"""

import numpy as np
import matplotlib.pyplot as plt

# 统一画图格式
plt.rcParams['font.size'] = 14
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'

# 读取数据文件
data = np.genfromtxt("lcurve.out", names=True)

# 绘制每个损失函数的损失值随训练步数的变化曲线
for name in data.dtype.names[1:-1]:
    plt.plot(data['step'], data[name], label=name)

# 设置坐标轴标签
plt.xlabel('Step')
plt.ylabel('Loss')

# 设置坐标轴的缩放比例
plt.xscale('symlog')
plt.yscale('log')

# 添加网格线和图例
plt.grid()
plt.legend()

# 保存图片为 PNG 文件
plt.savefig('./loss.png', bbox_inches='tight')
