# 脚本名称: process_force_data.py
# 运行方法: 在命令行或 IDE 中执行该脚本，确保提供正确的输入文件路径
# 脚本功能: 从给定的 .xyz 文件中读取结构数据，并根据力差异将其分类保存为两个不同的输出文件

from ase.io import read, write
import numpy as np
import os

def deal_data(file_name, force_data, fout='find_force_differ_10.xyz', fres='reserve_force_differ_10.xyz'):
    """
    处理输入的结构文件和力数据，根据力差异进行分类。

    参数:
        file_name (str): 输入的 .xyz 文件名，包含分子结构信息。
        force_data (str): 力数据文件名，每行对应一个结构的力信息。
        fout (str): 输出文件名，保存力差异超过阈值的结构。
        fres (str): 输出文件名，保存力差异在阈值内的结构。
    """

    # 读取结构文件和力数据
    test_frames = read(file_name, ":")
    data_force = np.loadtxt(force_data)

    # 检查力数据长度是否与原子总数匹配
    num_force_data = len(data_force)

    # 如果输出文件已存在，则先删除它们
    if os.path.exists(fout):
        os.remove(fout)

    if os.path.exists(fres):
        os.remove(fres)

    # 统计帧数和原子总数
    num_frames = 0
    total_atoms = 0

    # 遍历每个帧并累加原子数量
    for frame in test_frames:
        num_frames += 1
        total_atoms += frame.get_global_number_of_atoms()

    print(f"Number of frames: {num_frames} in {file_name}")
    print(f"Total number of atoms: {total_atoms} in {file_name}")

    if num_force_data != total_atoms:
        raise ValueError(f"The {file_name} frames do not match the force data in {force_data}")

    fout_frames, fres_frames, atom_ids = 0, 0, 0

    # 逐帧处理数据并检测力差异
    for frame_index, atoms in enumerate(test_frames):
        # 打印当前处理进度
        print(f"Processing frame {frame_index + 1} of {num_frames}")

        # 获取当前帧的原子数量
        num_atoms = atoms.get_global_number_of_atoms()
        next_frame_info = atom_ids + num_atoms

        # 获取 DFT 力数据
        DFT_force_inxyz = atoms.arrays["forces"].reshape(-1, 1)
        DFT_force_intrain = data_force[atom_ids:next_frame_info, 3:6].reshape(-1, 1)

        # 检查 DFT 数据是否与结构文件中的力数据匹配
        if not np.allclose(DFT_force_inxyz[-1], DFT_force_intrain[-1], atol=1e-06):
            raise ValueError(f"The {file_name} force does not match the force in {force_data}")

        # 确定预测力是否超过 DFT 力
        output_flag = get_outlier(
            data_force[atom_ids:next_frame_info, 0:3].reshape(-1, 1),
            data_force[atom_ids:next_frame_info, 3:6].reshape(-1, 1)
        )  # ev/A

        # 根据力差异将数据写入不同的文件
        if output_flag:
            write(fout, atoms, format='extxyz', append='ab')
            print('********* Outputting {}-st to {} *********'.format(frame_index + 1, fout))
            fout_frames += 1
        else:
            write(fres, atoms, format='extxyz', append='ab')
            # print('********* Outputting {}-st to {} *********'.format(frame_index+1, fres))
            fres_frames += 1

        atom_ids += num_atoms

    # 打印处理结果
    print(f"\nWe have a total of {num_frames} frames")
    print(f"Here we found {fout_frames} frames as outliers, and {fres_frames} as inliers")

def get_outlier(predicted_data, dft_data, diff_tol=10):
    """
    检查力差异是否超过给定的阈值。

    参数:
        predicted_data (np.ndarray): 预测的力数据。
        dft_data (np.ndarray): DFT 计算得到的力数据。
        diff_tol (float): 力差异的阈值，默认值为 10。

    返回:
        output_flag (bool): 如果力差异超过阈值，返回 True；否则返回 False。
    """
    force_outlier_atoms = 0
    output_flag = False

    # 遍历每个原子检查力差异
    for i in range(len(dft_data)):
        if abs(predicted_data[i] - dft_data[i]) > diff_tol:
            force_outlier_atoms += 1

    if force_outlier_atoms > 0:
        output_flag = True
        print(f"There are {force_outlier_atoms} outlier atoms in frame")

    return output_flag


if __name__ == "__main__":
    # 指定输入数据文件
    file_name = "train.xyz"
    force_data = "force_train.out"
    
    # 调用数据处理函数
    deal_data(file_name, force_data)
