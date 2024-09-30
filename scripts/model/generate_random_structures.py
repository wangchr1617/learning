# 文件名: generate_random_structures.py
# 运行方法: python generate_random_structures.py
# 功能描述: 从给定的结构文件生成随机扰动的原子结构，并将其保存到指定目录中。

import os
import numpy as np
from ase.build import sort
from ase.io import read, write
from RAG import random_atoms_gen as rag

def rag_iter(structure, idx, num_spec_dict, fix_ind_dict, iter_num=20, random_radi=0.01, cutoff_frac=0.9, strain_max=1.03, strain_min=0.97, vacuum_max=0.0):
    """
    使用RAG生成随机扰动结构并保存到文件。

    参数:
        structure (ase.Atoms): 输入结构。
        idx (str): 文件索引。
        num_spec_dict (dict): 原子数量的字典，格式如{'元素符号': 数量}。
        fix_ind_dict (dict): 固定原子索引的字典，格式如{'元素符号': [索引]}。
        iter_num (int, 可选): 迭代次数。默认为20。
        random_radi (float, 可选): 随机半径。默认为0.01。
        cutoff_frac (float, 可选): 截断分数。默认为0.9。
        strain_max (float, 可选): 最大应变。默认为1.03。
        strain_min (float, 可选): 最小应变。默认为0.97。
        vacuum_max (float, 可选): 最大真空。默认为0.0。
    """
    for i in range(iter_num):
        # 生成随机应变比例和真空
        strain_ratio = np.random.rand(3) * (strain_max - strain_min) + strain_min
        vacuum = [0., 0., np.random.rand() * vacuum_max]

        # 使用RAG生成随机结构
        atoms = rag(structure,
                    num_spec_dict=num_spec_dict,
                    fix_ind_dict=fix_ind_dict,
                    cutoff_frac=cutoff_frac,
                    random_radi=random_radi,
                    strain_ratio=strain_ratio,
                    vacuum=vacuum,
                    log=True
        )

        # 保存生成的结构
        output_file = f'./rag_structures/{idx}_{random_radi:.2f}-{i:02d}.vasp'
        write(output_file, sort(atoms))

def main():
    """
    主函数，生成指定结构的随机扰动。
    """
    # 确保输出目录存在
    rag_path = './rag_structures'
    if not os.path.exists(rag_path):
        os.makedirs(rag_path)

    # 定义骨架路径和参数
    backbone_path = './backbone/'
    radi_list = [0.08, 0.16, 0.32, 0.64]
    idx = '6'
    backbone = os.path.join(backbone_path, f'{idx}.vasp')

    # 读取输入结构
    structure = read(backbone)

    # 定义原子数量和固定索引
    num_spec_dict = {'Ge': 30, 'Te': 32}
    fix_ind = np.arange(30, 62)

    # 迭代生成不同随机半径的结构
    for random_radi in radi_list:
        fix_ind_dict = {'Te': fix_ind}
        rag_iter(structure, idx, num_spec_dict, fix_ind_dict, iter_num=20, random_radi=random_radi)

if __name__ == '__main__':
    main()
