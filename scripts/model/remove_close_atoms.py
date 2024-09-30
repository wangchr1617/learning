# 文件名: remove_close_atoms.py
# 运行方法: python remove_close_atoms.py
# 功能描述: 从XYZ文件中读取原子结构，删除距离小于给定阈值的原子对，并保存修改后的结构

from ase import Atoms
from ase.io import read, write
from scipy.spatial import cKDTree
import numpy as np

def remove_close_atoms(input_filename, output_filename, threshold):
    """
    从输入的XYZ文件中读取原子结构，删除距离小于给定阈值的原子对，并将修改后的结构保存到输出文件。

    参数:
    input_filename (str): 输入的XYZ文件名。
    output_filename (str): 输出的XYZ文件名。
    threshold (float): 判断原子对是否太近的距离阈值。
    """
    # 读取输入文件中的原子结构
    atoms = read(input_filename)
    # 获取原子的位置
    positions = atoms.get_positions()

    # 创建KD树以加速最近邻搜索
    tree = cKDTree(positions)

    # 查找距离小于阈值的原子对
    close_pairs = tree.query_pairs(threshold)

    # 初始化要删除的原子集合
    atoms_to_remove = set()

    # 遍历所有近距离原子对
    for i, j in close_pairs:
        # 选择要删除的原子
        atom_to_remove = min(i, j)
        atoms_to_remove.add(atom_to_remove)

    # 构建新的原子集合，不包含待删除的原子
    new_atoms = Atoms([atom for idx, atom in enumerate(atoms) if idx not in atoms_to_remove])

    # 将新的原子结构写入输出文件
    write(output_filename, new_atoms)

    print(f"Modified file saved as {output_filename}")

if __name__ == '__main__':
    # 输入文件名
    input_filename = "poly.xyz"
    # 输出文件名
    output_filename = "modified_poly.xyz"
    # 距离阈值
    threshold = 1.5

    # 调用函数处理原子结构
    remove_close_atoms(input_filename, output_filename, threshold)
