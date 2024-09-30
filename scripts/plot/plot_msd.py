# 文件名: plot_msd.py
# 运行方法: python plot_msd.py
# 功能描述: 读取均方位移（MSD）数据，并绘制MSD随时间变化的曲线图。

import numpy as np
import matplotlib.pyplot as plt

# 设置matplotlib全局参数
plt.rcParams['font.size'] = 10
plt.rcParams["axes.unicode_minus"] = False
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'

def load_msd_data(file_path):
    """
    加载MSD数据。

    参数:
        file_path (str): MSD数据文件路径。

    返回:
        tuple: 包含时间和MSD数据的元组。
    """
    data = np.loadtxt(file_path)

    correlation_time = data[:, 0]
    msd_x = data[:, 1]
    msd_y = data[:, 2]
    msd_z = data[:, 3]

    return correlation_time, msd_x, msd_y, msd_z

def plot_msd(correlation_time, msd_x, msd_y, msd_z):
    """
    绘制MSD随时间变化的曲线图。

    参数:
        correlation_time (array): 相关时间数组。
        msd_x (array): x方向的MSD数据。
        msd_y (array): y方向的MSD数据。
        msd_z (array): z方向的MSD数据。
    """
    plt.figure(figsize=(10, 6))
    plt.plot(correlation_time, msd_x, label='MSD in x direction')
    plt.plot(correlation_time, msd_y, label='MSD in y direction')
    plt.plot(correlation_time, msd_z, label='MSD in z direction')

    plt.xlabel('Correlation Time (ps)')
    plt.ylabel('MSD (Å²)')
    plt.title('MSD vs Correlation Time')
    plt.legend()

    # 设置x轴和y轴的范围
    plt.xlim(0, correlation_time[-1])
    plt.ylim(0, max(max(msd_x), max(msd_y), max(msd_z)) * 1.1)

    plt.tight_layout()
    plt.savefig("msd.png", bbox_inches='tight')
    print("MSD plot has been saved as msd.png.")

def main():
    """
    主函数，加载数据并调用绘图函数。
    """
    file_path = 'msd.out'
    correlation_time, msd_x, msd_y, msd_z = load_msd_data(file_path)
    plot_msd(correlation_time, msd_x, msd_y, msd_z)

if __name__ == '__main__':
    main()
