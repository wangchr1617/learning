# 文件名: plot_EV_nep.py
# 运行方法: python plot_EV_nep.py ev.txt nep.out
# 功能描述: 读取体积和DFT能量数据以及NEP能量数据，绘制能量-体积(EV)曲线，并标出能量最低点。

import matplotlib.pyplot as plt
import numpy as np
import sys

# 设置matplotlib全局参数
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'

def load_data(ev_file, nep_file):
    """
    从指定文件加载体积、DFT能量和NEP能量数据。

    参数:
        ev_file (str): 包含体积和DFT能量数据的文件路径。
        nep_file (str): 包含NEP能量数据的文件路径。

    返回:
        tuple: 包含体积、DFT能量和NEP能量的元组。
    """
    data = np.loadtxt(ev_file, delimiter=',', dtype=float)
    nep = np.loadtxt(nep_file, delimiter=' ', dtype=float)
    vol = data[:, 0]  # 体积数据
    dft = data[:, 1]  # DFT能量数据
    return vol, dft, nep

def find_min_energy_point(vol, dft):
    """
    找到能量最低点对应的体积和能量值。

    参数:
        vol (np.ndarray): 体积数组。
        dft (np.ndarray): DFT能量数组。

    返回:
        tuple: 能量最低点的体积和能量值。
    """
    min_index = np.argmin(dft)
    x = vol[min_index]
    y = dft[min_index]
    return x, y

def plot_ev_curve(vol, dft, nep, min_point):
    """
    绘制能量-体积(EV)曲线，并标出能量最低点。

    参数:
        vol (np.ndarray): 体积数组。
        dft (np.ndarray): DFT能量数组。
        nep (np.ndarray): NEP能量数组。
        min_point (tuple): 能量最低点的体积和能量值。
    """
    x, y = min_point
    point_label = f"({x:.2f}, {y:.2f})"

    plt.figure(figsize=(5, 3))
    plt.plot(vol, nep, linestyle='-.', color='orange', alpha=0.7, label="NEP")
    plt.plot(vol, dft, linestyle='--', alpha=0.7, label="DFT")
    plt.scatter(vol, dft, alpha=0.7)
    plt.scatter(x, y, color='red')  # 标记最低点
    plt.text(x, y + 0.1, point_label, verticalalignment="baseline", horizontalalignment="center")
    plt.xlabel(r'Volume / $Å^3$')
    plt.ylabel('Energy per atom / eV')
    plt.legend(loc="upper right")
    plt.tight_layout()  # 自动调整布局
    plt.savefig('./ev.png', bbox_inches='tight')
    print("EV plot has been saved as ev.png.")

def main():
    """
    主函数，解析命令行参数并绘制EV曲线。
    """
    if len(sys.argv) != 3:
        print("Usage: python EV_nep.py <ev_file> <nep_file>")
        sys.exit(1)

    # 加载数据
    vol, dft, nep = load_data(sys.argv[1], sys.argv[2])

    # 找到能量最低点
    min_point = find_min_energy_point(vol, dft)

    # 绘制EV曲线
    plot_ev_curve(vol, dft, nep, min_point)

if __name__ == '__main__':
    main()
