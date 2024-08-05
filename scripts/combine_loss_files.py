# 文件名: combine_loss_files.py
# 运行方法: python combine_loss_files.py loss1.out loss2.out
# 功能描述: 合并两个损失数据文件，并进行时间轴的平移。

import numpy as np
import sys

def load_loss_file(filename):
    """
    加载损失数据文件。

    参数:
        filename (str): 损失数据文件路径。

    返回:
        np.ndarray: 损失数据的数组。
    """
    try:
        loss_data = np.loadtxt(filename)
    except IOError:
        print(f"Error: Could not read file {filename}")
        sys.exit(1)
    except ValueError:
        print("Error: File format is not correct")
        sys.exit(1)
    return loss_data

def combine_loss_files(file1, file2):
    """
    合并两个损失数据文件，并进行时间轴的平移。

    参数:
        file1 (str): 第一个损失数据文件路径。
        file2 (str): 第二个损失数据文件路径。

    返回:
        np.ndarray: 合并后的损失数据。
    """
    # 加载损失数据
    loss1 = load_loss_file(file1)
    loss2 = load_loss_file(file2)

    # 时间轴平移
    loss2[:, 0] += loss1[-1, 0]

    # 合并数据
    combined_loss = np.concatenate((loss1, loss2), axis=0)

    return combined_loss

def save_combined_loss(combined_loss, output_file='loss_all.out'):
    """
    保存合并后的损失数据。

    参数:
        combined_loss (np.ndarray): 合并后的损失数据。
        output_file (str): 输出文件路径。
    """
    np.savetxt(output_file, combined_loss)

def main():
    """
    主函数，合并损失数据并保存结果。
    """
    if len(sys.argv) != 3:
        print("Usage: python combine_loss_files.py <loss1.out> <loss2.out>")
        sys.exit(1)

    # 合并损失数据
    combined_loss = combine_loss_files(sys.argv[1], sys.argv[2])

    # 保存合并后的损失数据
    save_combined_loss(combined_loss)

if __name__ == '__main__':
    main()
