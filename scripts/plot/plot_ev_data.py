"""
文件名: plot_ev_data.py
运行方法: python plot_ev_data.py
功能描述: 读取ev.txt文件中的数据并绘制能量-体积曲线。
"""

import matplotlib.pyplot as plt

# 设置画图格式
plt.rcParams['font.size'] = 14
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'

data_dict = {}

# 读取数据文件
with open('ev.txt', 'r') as file:
    label = None 
    for line in file:
        try:
            # 读取体积和能量数据
            vol, ene = map(float, line.split(','))
            if label is not None:
                data_dict[label]['volumes'].append(vol)
                data_dict[label]['energies'].append(ene)
        except ValueError:
            # 读取新的标签
            label = line.strip()  
            data_dict[label] = {'volumes': [], 'energies': []} 

# 创建绘图
plt.figure(figsize=(10, 6), dpi=300)
for label, data in data_dict.items():
    plt.plot(data['volumes'], data['energies'], linestyle='--', alpha=0.7, label=label)

# 设置坐标轴标签和图例
plt.xlabel(r'Volume / $Å^3$')
plt.ylabel('Energy per atom/eV')
plt.legend(loc="best")

# 保存图像
plt.savefig('./ev.png', bbox_inches='tight')
