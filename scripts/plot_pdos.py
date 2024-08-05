# 文件名: plot_pdos.py
# 运行方法: python plot_pdos.py PDOS_Sb.dat
# 功能描述: 读取PDOS数据文件，并绘制投影态密度（PDOS）图。

import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import sys

# 设置plt的参数
plt.rcParams['font.size'] = 14
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'

def plot_pdos(file_path):
    """
    读取PDOS数据文件并绘制投影态密度（PDOS）图。

    参数:
    file_path (str): PDOS数据文件的路径。
    """
    # 设置输出图像的文件名
    target = os.path.splitext(file_path)[0] + '.png'

    # 读取文件并解析数据
    with open(file_path, 'r') as f:
        data = []
        for line in f.readlines():
            # 去除换行符并按空格分隔数据项
            line_data = [item for item in line.replace("\n", "").split(" ") if item != ""]
            data.append(line_data)

    # 将数据转换为DataFrame
    df = pd.DataFrame(data)

    # 提取x轴数据并转换为数值类型
    x = pd.to_numeric(df[0][1:].reset_index(drop=True), errors='coerce')

    # 初始化y轴最大值列表
    y_max = []

    # 绘制每一列的数据
    for i in range(4):
        # 提取y轴数据并转换为数值类型
        y = pd.to_numeric(df[i + 1][1:].reset_index(drop=True), errors='coerce')
        # 绘制曲线
        plt.plot(x, y, alpha=1, linewidth=1, label=df[i + 1][0])
        # 记录y轴的最大值
        y_max.append(y.max())

    # 设置x轴和y轴标签
    plt.xlabel('Energy (eV)')
    plt.ylabel('DOS (states/eV)')

    # 设置x轴和y轴的范围
    plt.xlim(x.min(), x.max())
    plt.ylim(0, max(y_max) * 1.1)

    # 显示图例
    plt.legend(loc='best')

    # 保存图像
    plt.savefig(target, bbox_inches='tight')

    # 输出保存文件名
    print(f"PDOS图像已保存为: {target}")

if __name__ == '__main__':
    # 获取命令行参数中的文件路径
    path = sys.argv[1]

    # 调用绘图函数
    plot_pdos(path)
