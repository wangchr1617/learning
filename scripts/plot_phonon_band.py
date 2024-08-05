# 文件名: plot_phonon_band.py
# 运行方法: python plot_phonon_band.py band.dat
# 功能描述: 读取声子频率的band.dat文件并绘制声子能带图

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys

# 设置plt的参数
plt.rcParams['font.size'] = 14
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'

def main():
    # 读取命令行参数中的band.dat文件
    input_file = sys.argv[1]

    # 打开文件并读取所有行
    with open(input_file) as f:
        lines = f.readlines()

    # 读取文件头部的注释
    comment = lines.pop(0)
    
    # 读取k点的标签
    klabels = [eval(i) for i in lines.pop(0).strip().split(" ")[3:]]
    
    # 初始化数据容器
    data = []
    band = []
    
    # 逐行读取数据
    for line in lines:
        row = line.strip().split(" ")
        
        # 检查空行，用于分隔不同的数据块
        if "" in row:
            data.append(band)
            band = []
            continue
        
        # 读取每一行的数据
        band.append([eval(row[0]), eval(row[1])])
    
    # 创建图形和坐标轴对象
    fig, ax = plt.subplots(figsize=(4.2, 3), dpi=140)

    # 绘制声子能带图
    for band_data in data:
        df = pd.DataFrame(band_data)
        
        # 提取x和y坐标
        x = df[0]
        y = df[1]
        
        # 绘制每个能带
        ax.plot(x, y, c='r', lw=0.8)

    # 设置x轴的范围和标签
    ax.set_xlim(min(klabels), max(klabels))
    ax.set_ylabel('Frequency (THz)')
    ax.set_xticks(klabels)

    # 设置k点的标签（可以根据需要取消注释并自定义标签）
    # BAND_LABELS = [r'$\Gamma$', 'T', r'H$_2$', r'H$_0$', 'L', r'$\Gamma$', r'S$_0$', r'S$_2$', 'F', r'$\Gamma$']
    # ax.set_xticklabels(BAND_LABELS)

    # 绘制垂直虚线和水平虚线
    for x in klabels:
        ax.axvline(x, c='0.5', ls="--", lw=0.8)
    ax.axhline(y=0.0, c='0.5', ls="--", lw=0.8)

    # 调整布局以适应图形
    plt.tight_layout()

    # 保存图形到文件
    plt.savefig('./phon.png', bbox_inches='tight')

    # 输出完成提示
    print("声子能带图已保存为phon.png")

if __name__ == "__main__":
    main()
