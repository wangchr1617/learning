"""
用法: python rand_vacancy.py

功能描述:
该脚本从VASP的POSCAR文件中随机删除指定数量的特定类型的原子，生成具有空位的结构，并保存为新的VASP文件。

参数说明:
- `filename`: 输入的VASP文件名，默认读取当前目录下的 `POSCAR` 文件。
- `A`: 需要删除的原子类型。
- `n`: 删除的原子数量。

结果:
- 生成具有指定数量空位的VASP文件，文件名格式为 `vac_[n].vasp`。
"""

from ase.io import read, write
from ase.build import sort
import random

def region(structure, X):
    """
    获取指定类型原子在结构中的起始和结束索引。

    参数:
    - structure: ASE Atoms 对象。
    - X: 目标原子类型的符号。

    返回:
    - begin: 起始索引。
    - end: 结束索引。
    """
    AtomList = structure.get_chemical_symbols()
    begin = AtomList.index(X)
    for i in range(begin, len(AtomList)):
        if AtomList[i] != X:
            end = i
            break
    else:
        end = len(AtomList)
    return begin, end

def rand_vacancy(filename, A, n):
    """
    在给定结构中随机删除指定数量的原子。

    参数:
    - filename: 输入的VASP文件名。
    - A: 需要删除的原子类型。
    - n: 删除的原子数量。

    输出:
    - 生成新的VASP文件，文件名格式为 `vac_[n].vasp`。
    """
    structure = sort(read(filename))  # 读取并排序结构
    begin, end = region(structure, A)  # 获取目标原子的索引范围

    if end - begin < n:
        print("Error: Not enough atoms of type '{}' to remove {}.".format(A, n))
        return

    # 随机选择要删除的原子索引
    Alist = random.sample(range(begin, end), n)
    print(f"Removing atoms at indices: {Alist}")

    # 删除选择的原子并生成新结构
    del structure[Alist]
    output_filename = f"vac_{n}.vasp"
    write(output_filename, sort(structure))
    print(f"Vacancy structure saved as: {output_filename}")

if __name__ == "__main__":
    rand_vacancy('./POSCAR', 'Ge', 150)
