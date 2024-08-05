# 文件名: exchangeABatoms.py
# 运行方法: python exchangeABatoms.py POSCAR Ge Te 1 #001
# 功能描述: 在给定的POSCAR文件中，随机交换两个元素的原子位置，并将结果保存到新的文件中。

from ase.io import read, write
from ase.build import sort
import os
import random
import sys

def region(structure, element):
    """
    确定指定元素在结构中的区域（开始和结束索引）。

    参数:
        structure (Atoms): ASE中的原子结构对象。
        element (str): 目标元素符号。

    返回:
        tuple: 元素在结构中的开始和结束索引。
    """
    atom_list = structure.get_chemical_symbols()
    begin = atom_list.index(element)
    for i in range(begin, len(atom_list)):
        if atom_list[i] != element:
            end = i
            break
    else:
        end = len(atom_list)
    return begin, end

def rand_exchange(filename, element_a, element_b, n=1, idx=None):
    """
    随机交换两个元素的原子位置。

    参数:
        filename (str): 输入POSCAR文件名。
        element_a (str): 元素A的符号。
        element_b (str): 元素B的符号。
        n (int): 交换的原子数。
        idx (str, optional): 输出文件名中的索引。
    """
    # 读取并排序结构
    structure = sort(read(filename))
    name = os.path.splitext(filename)[0]

    # 确定元素A和B的位置范围
    a_begin, a_end = region(structure, element_a)
    b_begin, b_end = region(structure, element_b)

    # 检查是否有足够的原子进行交换
    if a_end - a_begin < n or b_end - b_begin < n:
        print("Error!!!")
        print("The atom number of element is less than the number you want to exchange!")
        return
    
    # 创建输出目录
    if not os.path.exists('./exchange'):
        os.makedirs('./exchange')

    # 随机选择要交换的原子索引
    a_list = random.sample(range(a_begin, a_end), n)
    b_list = random.sample(range(b_begin, b_end), n)

    # 执行交换
    for i in a_list:
        structure[i].symbol = element_b
    for i in b_list:
        structure[i].symbol = element_a

    # 保存交换后的结构
    if idx is None:
        idx = n
    output_filename = f"./exchange/{name}_{idx}.vasp"
    write(output_filename, sort(structure))
    print(f"Exchange completed. Output saved to {output_filename}")

def main():
    """
    主函数，解析命令行参数并调用交换函数。
    """
    if len(sys.argv) < 5:
        print("Usage: python exchangeABatoms.py POSCAR Ge Te 1 #001")
        sys.exit(1)

    # 获取命令行参数
    structure = sys.argv[1]
    element_a = sys.argv[2]
    element_b = sys.argv[3]
    n = int(sys.argv[4])

    # 检查是否提供了索引参数
    idx = None
    if len(sys.argv) > 5:
        idx = sys.argv[5]

    # 调用交换函数
    rand_exchange(structure, element_a, element_b, n, idx)

if __name__ == "__main__":
    main()
