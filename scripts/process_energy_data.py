'''
脚本名称: process_energy_data.py
运行方法: 在命令行或 IDE 中执行该脚本，确保提供正确的输入文件路径
脚本功能: 从给定的 .xyz 文件中读取结构数据，并根据能量差异将其分类保存为两个不同的输出文件
'''

from ase.io import read, write
import numpy as np
import os

def deal_data(file_name, energy_data, fout='find_energy_differ_0.5eV.xyz', fres='reserve_energy_differ_0.5eV.xyz'):
    """
    处理输入的结构文件和能量数据，根据能量差异进行分类。

    参数:
        file_name (str): 输入的 .xyz 文件名，包含分子结构信息。
        energy_data (str): 能量数据文件名，每行对应一个结构的能量。
        fout (str): 输出文件名，保存能量差异超过阈值的结构。
        fres (str): 输出文件名，保存能量差异在阈值内的结构。
    """
    
    # 读取结构文件和能量数据
    test_frames = read(file_name, ":")
    data_energy = np.loadtxt(energy_data)

    # 检查能量数据长度
    num_energy_data = len(data_energy)

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

    fout_frames, fres_frames = 0, 0

    # 逐帧处理数据并检测能量差异
    for frame_index, atoms in enumerate(test_frames):
        # 打印当前处理进度
        print(f"Processing frame {frame_index + 1} of {num_frames}")
        
        # 获取当前帧的原子数量
        num_atoms = atoms.get_global_number_of_atoms()

        # 计算每个原子的能量
        energy = atoms.info["energy"] / num_atoms
        
        # 从能量数据中获取相应的 DFT 能量
        DFT_energy_intrain = data_energy[frame_index, 1]

        # 检查 DFT 数据是否与结构文件中的能量匹配
        if not np.allclose(energy, DFT_energy_intrain, atol=1e-06):
            raise ValueError(f"The {file_name} energy does not match the energy in {energy_data}")

        # 确定预测能量是否超过 DFT 能量
        output_flag = get_outlier(data_energy[frame_index, 0], data_energy[frame_index, 1]) # eV/atom

        # 根据能量差异将数据写入不同的文件
        if output_flag:
            write(fout, atoms, format='extxyz', append='ab')
            print('********* Outputting {}-st to {} *********'.format(frame_index + 1, fout))
            fout_frames += 1
        else:
            write(fres, atoms, format='extxyz', append='ab')
            # print('********* Outputting {}-st to {} *********'.format(frame_index+1, fres))
            fres_frames += 1

    # 打印处理结果
    print(f"We have a total of {num_frames} frames")
    print(f"\nWe found {fout_frames} frames as outliers, and {fres_frames} as inliers")

def get_outlier(predicted_data, dft_data, diff_tol=0.5):
    """
    检查能量差异是否超过给定的阈值。

    参数:
        predicted_data (float): 预测的能量数据。
        dft_data (float): DFT 计算得到的能量数据。
        diff_tol (float): 能量差异的阈值，默认值为 0.5。

    返回:
        output_flag (bool): 如果能量差异超过阈值，返回 True；否则返回 False。
    """
    output_flag = False
    
    if abs(predicted_data - dft_data) > diff_tol:
        output_flag = True
        print(f"Energy outlier detected for this frame")

    return output_flag


if __name__ == "__main__":
    # 指定输入数据文件
    file_name = "train.xyz"
    energy_data = "energy_train.out"
    
    # 调用数据处理函数
    deal_data(file_name, energy_data)
