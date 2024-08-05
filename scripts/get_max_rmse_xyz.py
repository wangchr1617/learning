# 文件名: get_max_rmse_xyz.py
# 运行方法:
#     python get_max_rmse_xyz.py train.xyz force_train.out 13
#     python get_max_rmse_xyz.py train.xyz virial_train.out 13 
#     python get_max_rmse_xyz.py train.xyz energy_train.out 13 
# 功能描述: 获取force_train.out、virial_train.out或energy_train.out文件中误差最大的点，
# 并找到这些点所属的训练集ID并输出。
# 注意: 输出的轨迹数小于或等于误差点的数量。

import numpy as np
import sys

def get_rmse_ids(nmax, file_force_loss):
    """
    获取误差最大的点的RMSE和ID。

    参数:
        nmax (int): 要获取的最大误差点的数量。
        file_force_loss (str): 包含力误差数据的文件路径。

    返回:
        tuple: 包含最大RMSE和对应ID的元组。
    """
    frmse = np.loadtxt(file_force_loss)
    # 判断数据格式并计算RMSE
    if frmse.shape[1] == 6:
        rmse = np.sum(np.abs(frmse[:, 0:3] - frmse[:, 3:6]), axis=1)
    else:
        rmse = np.abs(frmse[:, 0] - frmse[:, 1])
    rmse_max_ids = np.argsort(-rmse)

    return rmse[rmse_max_ids[:nmax]], rmse_max_ids[:nmax]

def get_frame_lines(train_xyz):
    """
    获取训练集文件中的每个帧的起始行和原子数。

    参数:
        train_xyz (str): 训练集文件路径。

    返回:
        tuple: 包含帧起始行号列表和原子数列表的元组。
    """
    num_lines, num_atoms = [], []
    with open(train_xyz, "r") as fi:
        flines = fi.readlines()
        for i, line in enumerate(flines):
            if "energy" in line or "Energy" in line:
                num_lines.append(i - 1)
                num_atoms.append(int(flines[i - 1]))

    return num_lines, num_atoms

def print_max_xyz(frame_list, num_lines, train_xyz, fout="find_out.xyz"):
    """
    将最大误差点所属的帧写入新的XYZ文件。

    参数:
        frame_list (list): 帧ID列表。
        num_lines (list): 帧起始行号列表。
        train_xyz (str): 训练集文件路径。
        fout (str): 输出文件名。
    """
    fout_str = ""
    with open(train_xyz, "r") as fi, open(fout, 'w') as fo:
        flines = fi.readlines()
        for frame_id in frame_list:
            flsta = num_lines[frame_id]
            if frame_id == len(num_lines) - 1:
                flend = len(flines)
            else:
                flend = num_lines[frame_id + 1]
            fout_str += "".join(flines[flsta:flend])

        fo.write(fout_str)

def main():
    """
    主函数，解析命令行参数并执行查找最大误差点的逻辑。
    """
    if len(sys.argv) != 4:
        print("Usage: python get_max_rmse_xyz.py <train_xyz> <force_loss_file> <nmax>")
        sys.exit(1)

    fxyz = sys.argv[1]
    floss = sys.argv[2]
    nmax = int(sys.argv[3])

    # 获取最大误差点的RMSE和ID
    rmse_max, rmse_ids = get_rmse_ids(nmax, floss)
    num_lines, num_atoms = get_frame_lines(fxyz)
    sum_atoms = np.cumsum(num_atoms)

    frame_list = []
    for i in rmse_ids:
        nframe = np.searchsorted(sum_atoms, i)
        if nframe not in frame_list:
            frame_list.append(nframe)

    print(f"The largest RMSE with {nmax} atoms are located in frame {frame_list}")
    print_max_xyz(frame_list, num_lines, fxyz, fout='find_out.xyz')

if __name__ == '__main__':
    main()
