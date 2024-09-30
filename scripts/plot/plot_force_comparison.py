# 文件名: plot_force_comparison.py
# 运行方法: python plot_force_comparison.py
# 功能描述: 读取计算的力数据和损失数据，绘制DFT和NEP力的比较图，并计算训练集的RMSE。

import matplotlib.pyplot as plt
import numpy as np

# 设置matplotlib全局参数
plt.rcParams['font.size'] = 14
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'

def load_data(force_file, loss_file):
    """
    加载力数据和损失数据。

    参数:
        force_file (str): 力数据文件路径。
        loss_file (str): 损失数据文件路径。

    返回:
        tuple: 包含力数据和损失数据的numpy数组。
    """
    try:
        force_data = np.loadtxt(force_file)
        loss_data = np.loadtxt(loss_file)
    except IOError as e:
        print(f"Error: Could not read file: {e}")
        sys.exit(1)
    except ValueError as e:
        print(f"Error: File format is not correct: {e}")
        sys.exit(1)
    return force_data, loss_data

def plot_force_comparison(force, loss, output_file='force.png'):
    """
    绘制DFT和NEP力的比较图。

    参数:
        force (np.ndarray): 力数据的数组。
        loss (np.ndarray): 损失数据的数组。
        output_file (str): 输出图像文件路径。
    """
    # 提取NEP和DFT力数据
    nep = force[:, 0:3]
    dft = force[:, 3:6]

    # 确定图像坐标范围
    lim0, lim1 = int(dft.min()) - 10, int(dft.max()) + 10
    distance = lim1 - lim0

    # 创建图形并绘制数据点
    plt.figure(figsize=(5, 4))
    plt.plot(dft, nep, 'o', markersize=3)
    plt.plot(np.linspace(lim0, lim1), np.linspace(lim0, lim1), '-', c='#85618E')
    plt.xlim(lim0, lim1)
    plt.ylim(lim0, lim1)
    plt.xlabel('DFT force (eV/Å)')
    plt.ylabel('NEP force (eV/Å)')

    # 在图中添加RMSE文本
    string = 'train RMSE = {:.1f} meV/Å'.format(loss[-1, 5] * 1000)
    plt.text(lim1 - distance * 0.05, lim0 + distance * 0.05, string, horizontalalignment='right', verticalalignment='baseline')

    # 添加图例
    plt.legend(['x direction', 'y direction', 'z direction'])

    # 调整布局并保存图像
    plt.tight_layout()
    plt.savefig(output_file, bbox_inches='tight')

def main():
    """
    主函数，加载数据并绘制图形。
    """
    # 加载力数据和损失数据
    force, loss = load_data('./force_train.out', './loss.out')

    # 绘制力比较图
    plot_force_comparison(force, loss)

if __name__ == '__main__':
    main()
