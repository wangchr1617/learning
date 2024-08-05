# 文件名: plot_energy_comparison.py
# 运行方法: python plot_energy_comparison.py
# 功能描述: 读取计算的能量数据和损失数据，绘制DFT和NEP能量的比较图，并计算训练集的RMSE。

import matplotlib.pyplot as plt
import numpy as np

# 设置matplotlib全局参数
plt.rcParams['font.size'] = 14
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'

def load_data(energy_file, loss_file):
    """
    加载能量数据和损失数据。

    参数:
        energy_file (str): 能量数据文件路径。
        loss_file (str): 损失数据文件路径。

    返回:
        tuple: 包含能量数据和损失数据的numpy数组。
    """
    try:
        energy_data = np.loadtxt(energy_file)
        loss_data = np.loadtxt(loss_file)
    except IOError as e:
        print(f"Error: Could not read file: {e}")
        sys.exit(1)
    except ValueError as e:
        print(f"Error: File format is not correct: {e}")
        sys.exit(1)
    return energy_data, loss_data

def plot_energy_comparison(energy, loss, output_file='energy.png'):
    """
    绘制DFT和NEP能量的比较图。

    参数:
        energy (np.ndarray): 能量数据的数组。
        loss (np.ndarray): 损失数据的数组。
        output_file (str): 输出图像文件路径。
    """
    # 提取NEP和DFT能量数据
    nep = energy[:, 0]
    dft = energy[:, 1]

    # 确定图像坐标范围
    lim0, lim1 = int(dft.min()) - 1, int(dft.max()) + 1
    distance = lim1 - lim0

    # 创建图形并绘制数据点
    plt.figure(figsize=(5, 4))
    plt.plot(dft, nep, 'o', c='#E57B71', markersize=3, label='Energy Data')
    plt.plot(np.linspace(lim0, lim1), np.linspace(lim0, lim1), '-', c='#85618E', label='y=x Line')
    plt.xlim(lim0, lim1)
    plt.ylim(lim0, lim1)
    plt.xlabel('DFT energy (eV/atom)')
    plt.ylabel('NEP energy (eV/atom)')

    # 在图中添加RMSE文本
    string = 'train RMSE = {:.1f} meV/atom'.format(loss[-1, 4] * 1000)
    plt.text(lim1 - distance * 0.05, lim0 + distance * 0.05, string, horizontalalignment='right', verticalalignment='baseline')

    # 添加图例
    plt.legend()

    # 调整布局并保存图像
    plt.tight_layout()
    plt.savefig(output_file, bbox_inches='tight')

def main():
    """
    主函数，加载数据并绘制图形。
    """
    # 加载能量数据和损失数据
    energy, loss = load_data('./energy_train.out', './loss.out')

    # 绘制能量比较图
    plot_energy_comparison(energy, loss)

if __name__ == '__main__':
    main()
