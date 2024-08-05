# 文件名: plot_loss.py
# 运行方法: python plot_loss.py
# 功能描述: 读取损失数据文件，并绘制每一代的损失函数变化图。

import matplotlib.pyplot as plt
import numpy as np

# 设置matplotlib全局参数
plt.rcParams['font.size'] = 14
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'

def load_loss_data(file_path):
    """
    从文件中加载损失数据。

    参数:
        file_path (str): 损失数据文件的路径。

    返回:
        np.ndarray: 损失数据的数组。
    """
    try:
        return np.loadtxt(file_path)
    except IOError:
        print(f"Error: Could not read file {file_path}")
        sys.exit(1)
    except ValueError:
        print("Error: File format is not correct")
        sys.exit(1)

def plot_loss(loss, columns, output_file='loss.png'):
    """
    绘制损失函数变化图。

    参数:
        loss (np.ndarray): 损失数据的数组。
        columns (list): 损失数据的列标签。
        output_file (str): 输出图像文件路径。
    """
    plt.figure(figsize=(5, 4))
    # 绘制各个损失成分的变化曲线
    plt.plot(loss[:, 0], loss[:, 1:7])
    # 使用对数坐标轴
    plt.loglog()
    plt.xlim(100, loss[-1, 0])
    plt.xlabel('Generation')
    plt.ylabel('Loss')
    plt.legend(columns[1:7])
    plt.tight_layout()
    plt.savefig(output_file, bbox_inches='tight')

def main():
    """
    主函数，加载损失数据并绘制图形。
    """
    # 定义损失列标签
    columns = [
        'Generation', 'Total', 'L1-regularization', 'L2-regularization', 'Energy-train', 'Force-train',
        'Virial-train', 'Energy-test', 'Force-test', 'Virial-test'
    ]

    # 加载损失数据
    loss = load_loss_data('loss.out')

    # 绘制损失图
    plot_loss(loss, columns)

if __name__ == '__main__':
    main()
