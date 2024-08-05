# 文件名: plot_EV.py
# 运行方法: python plot_EV.py ev.txt
# 功能描述: 读取包含体积和能量数据的文本文件，绘制能量-体积(EV)曲线，并在图中标出能量最低点。

import matplotlib.pyplot as plt
import numpy as np
import sys

def load_data(file_path):
    """
    从指定文件加载体积和能量数据。

    参数:
        file_path (str): 包含体积和能量数据的文件路径。

    返回:
        tuple: 包含体积数组和能量数组的元组。
    """
    data = np.loadtxt(file_path, delimiter=',', dtype=float)
    vol = data[:, 0]  # 体积数据
    dft = data[:, 1]  # 能量数据
    return vol, dft

def find_min_energy_point(vol, dft):
    """
    找到能量最低点对应的体积和能量值。

    参数:
        vol (np.ndarray): 体积数组。
        dft (np.ndarray): 能量数组。

    返回:
        tuple: 能量最低点的体积和能量值。
    """
    min_index = np.argmin(dft)
    x = vol[min_index]
    y = dft[min_index]
    return x, y

def plot_ev_curve(vol, dft, min_point):
    """
    绘制能量-体积(EV)曲线，并标出能量最低点。

    参数:
        vol (np.ndarray): 体积数组。
        dft (np.ndarray): 能量数组。
        min_point (tuple): 能量最低点的体积和能量值。
    """
    x, y = min_point
    point_label = f"({x:.2f}, {y:.2f})"

    plt.figure(figsize=(5, 3))
    plt.plot(vol, dft, linestyle='--', alpha=0.7, label="DFT")
    plt.scatter(vol, dft, alpha=0.7)
    plt.scatter(x, y, c='r')  # 标记最低点
    plt.text(x, y + 0.05, point_label, verticalalignment="baseline", horizontalalignment="center")
    plt.xlabel(r'Volume / $Å^3$')
    plt.ylabel('Energy per atom / eV')
    plt.legend(loc="best")
    plt.tight_layout()  # 自动调整布局
    plt.savefig('./ev.png', bbox_inches='tight')
    print("EV plot has been saved as ev.png.")

def main():
    """
    主函数，解析命令行参数并绘制EV曲线。
    """
    if len(sys.argv) != 2:
        print("Usage: python EV_plot.py <ev_file>")
        sys.exit(1)

    # 加载数据
    vol, dft = load_data(sys.argv[1])

    # 找到能量最低点
    min_point = find_min_energy_point(vol, dft)

    # 绘制EV曲线
    plot_ev_curve(vol, dft, min_point)

if __name__ == '__main__':
    main()
