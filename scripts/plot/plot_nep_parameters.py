# 文件名: plot_nep_parameters.py
# 运行方法: python plot_nep_parameters.py
# 功能描述: 读取nep.restart文件中的参数，并绘制参数值随索引变化的图像。

import matplotlib.pyplot as plt
import numpy as np

# 设置matplotlib全局参数
plt.rcParams['font.size'] = 14
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'

def load_parameters(file_path):
    """
    从文件中加载参数数据。

    参数:
        file_path (str): 参数文件的路径。

    返回:
        np.ndarray: 参数数据的数组。
    """
    return np.loadtxt(file_path)

def plot_parameters(parameters):
    """
    绘制参数值随索引变化的图像。

    参数:
        parameters (np.ndarray): 参数数据数组。
    """
    x = np.arange(len(parameters))
    y = parameters[:, 0]

    plt.figure(figsize=(5, 4))
    plt.plot(x, y, label='Parameter 0')
    plt.xlim(min(x), max(x))
    plt.xlabel('Index')
    plt.ylabel('Parameter Value')
    plt.title('NEP Parameters Plot')
    plt.legend()
    plt.tight_layout()
    plt.savefig('./para.png', bbox_inches='tight')

def main():
    """
    主函数，加载参数数据并绘制图形。
    """
    parameters = load_parameters('./nep.restart')
    plot_parameters(parameters)

if __name__ == '__main__':
    main()
