# 文件名: plot_virial_comparison.py
# 运行方法: python plot_virial_comparison.py
# 功能描述: 加载数据并绘制DFT和NEP virial的比较图。

import matplotlib.pyplot as plt
import numpy as np

# 设置plt的参数
plt.rcParams['font.size'] = 14
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'

def load_data():
    """
    加载数据文件。

    返回:
        tuple: 包含virial和loss数据的numpy数组。
    """
    virial = np.loadtxt('./virial_train.out')
    loss = np.loadtxt('./loss.out')
    return virial, loss

def plot_virial_comparison(virial, loss):
    """
    绘制DFT和NEP virial的比较图。

    参数:
        virial (np.ndarray): virial数据。
        loss (np.ndarray): 损失数据。
    """
    # 提取NEP和DFT virial数据
    nep = virial[:, 0:6]
    dft = virial[:, 6:12]

    # 确定图像坐标范围
    lim0, lim1 = int(dft.min()) - 1, int(dft.max()) + 1
    distance = lim1 - lim0

    # 创建图形并绘制数据点
    plt.figure(figsize=(5, 4))
    plt.plot(dft, nep, 'o', markersize=3)
    plt.plot(np.linspace(lim0, lim1), np.linspace(lim0, lim1), '-', c='#85618E')
    plt.xlim(lim0, lim1)
    plt.ylim(lim0, lim1)
    plt.xlabel('DFT virial (eV/atom)')
    plt.ylabel('NEP virial (eV/atom)')

    # 在图中添加RMSE文本
    string = 'train RMSE = {:.1f} meV/Å'.format(loss[-1, 6] * 1000)
    plt.text(lim1 - distance * 0.05, lim0 + distance * 0.05, string, horizontalalignment='right', verticalalignment='baseline')

    # 添加图例
    plt.legend(['xx', 'yy', 'zz', 'xy', 'yz', 'zx'])

    # 调整布局并保存图像
    plt.tight_layout()
    plt.savefig('./virial.png', bbox_inches='tight')

def main():
    """
    主函数，加载数据并绘制图形。
    """
    virial, loss = load_data()
    plot_virial_comparison(virial, loss)

if __name__ == '__main__':
    main()
