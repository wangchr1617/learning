# 文件名: generate_random_perturbed_structures.py
# 运行方法: python generate_random_perturbed_structures.py
# 功能描述: 根据给定的POSCAR文件生成具有随机晶格扰动的结构文件

import os
import numpy as np
from pymatgen.io.vasp.inputs import Poscar
from pymatgen.core.lattice import Lattice

def generate_random_perturbed_structures(filename, tolerance=None, n_struct=None):
    """
    生成具有随机晶格扰动的结构文件。
    
    参数:
    filename (str): 原始结构文件的路径。
    tolerance (list): 改变晶格长度和角度的容差，是一个包含两个整数的列表。
                      默认值为[3, 3]。
    n_struct (int): 要生成的结构数量。默认值为10。
    """
    # 默认参数设置
    if tolerance is None:
        tolerance = [3, 3]

    if n_struct is None:
        n_struct = 10

    # 读取原始结构
    origin_struct = Poscar.from_file(filename, check_for_POTCAR=False).structure

    # 获取原始晶格参数
    old_latt = origin_struct.lattice.as_dict(verbosity=1)
    old_latt = list(old_latt.values())[4:10]

    # 初始化变量
    deform_struct = origin_struct
    delta, sigma = map(int, tolerance)
    Dlength = range(-delta, delta+1)
    Dangle = range(-sigma, sigma+1)
    lengths = old_latt[:3]
    angles = old_latt[3:]

    # 创建保存扰动结构的目录
    Pfile = 'perturb'
    if not os.path.exists(Pfile):
        os.makedirs(Pfile)

    # 生成具有随机扰动的结构
    for rnd in range(n_struct):
        lengths_loss = np.random.choice(Dlength, 3) / 100  # 随机选择长度扰动
        lengths_loss = np.array(lengths) * lengths_loss
        angle_loss = np.random.choice(Dangle, 3)  # 随机选择角度扰动
        loss = np.append(lengths_loss, angle_loss)
        Perturbed_lattice = np.array(old_latt) - loss
        
        # 使用新的晶格参数创建新晶格
        a, b, c, alpha, beta, gamma = Perturbed_lattice
        new_latt = Lattice.from_parameters(a, b, c, alpha, beta, gamma)
        deform_struct.lattice = new_latt
        
        # 将扰动结构写入文件
        out_poscar = Poscar(deform_struct, comment=f'{rnd}_Perturbed structure')
        out_poscar.write_file(f'{Pfile}/perturbed_{rnd}.vasp')

    print(f"{n_struct} structures with random perturbation in lattice were generated! Bye!")

if __name__ == '__main__':
    # 设置参数
    file = 'POSCAR'  # 输入文件名
    tolerance = [3, 3]  # 长度和角度扰动的容差
    size = 50  # 生成结构的数量
    
    # 调用函数生成结构
    generate_random_perturbed_structures(file, tolerance, size)
