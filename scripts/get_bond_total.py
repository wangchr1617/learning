# 文件名: get_bond_total.py
# 运行方法: python get_bond_total.py
# 功能描述: 从POSCAR文件中读取结构数据，计算所有可能的原子对之间的距离，并输出结果。

import os
from collections import defaultdict
from itertools import combinations

import numpy as np
from pymatgen.core.structure import Structure

def calculate_pair_distances(structure):
    """
    计算给定结构中所有可能的原子对之间的距离。

    参数:
        structure (Structure): pymatgen中的结构对象。

    返回:
        list: 包含所有原子对距离和对应的原子对名称的列表。
    """
    sites = structure.sites
    indices = [idx + 1 for idx, site in enumerate(sites)]  # 因为Python的索引编号是从0开始的

    # 获取所有可能的原子对和对应的索引对
    pairs = list(combinations(sites, r=2))
    idxs = list(combinations(indices, r=2))

    # 计算每个原子对的距离并存储
    pairs_dist = defaultdict(list)
    for pair, idx in zip(pairs, idxs):
        d = structure.lattice.get_all_distances(pair[0].frac_coords, pair[1].frac_coords)[0][0]
        pair_name = f"{pair[0].specie}{idx[0]}-{pair[1].specie}{idx[1]}"
        pairs_dist[np.round(d, decimals=3)].append(pair_name)

    # 按距离排序
    return sorted(pairs_dist.items())

def print_structure_info(structure):
    """
    打印结构信息，包括原子编号和种类。

    参数:
        structure (Structure): pymatgen中的结构对象。
    """
    print("Note: --------------------")
    print("{:<8} {}".format("Number", "Element"))
    for idx, site in enumerate(structure.sites):
        print("{:<8} {:<8}".format(idx + 1, site.specie))

def print_distance_info(pairs_dist):
    """
    打印所有原子对的距离信息。

    参数:
        pairs_dist (list): 包含所有原子对距离和对应的原子对名称的列表。
    """
    print("Note: --------------------")
    print("{:<8} {:<8} {}".format("Species", "Species", "Distance"))
    for d, pairs in pairs_dist:
        print(f"distance = {d}, Number = {len(pairs)}")
        print("    " + "   ".join(pairs))

def save_distance_info(pairs_dist, filename="distance.dat"):
    """
    将所有原子对的距离信息保存到文件。

    参数:
        pairs_dist (list): 包含所有原子对距离和对应的原子对名称的列表。
        filename (str): 保存的文件名。
    """
    with open(filename, "w") as f:
        for d, pairs in pairs_dist:
            f.write(f"distance = {d}, Number = {len(pairs)}\n")
            f.write("    " + "   ".join(pairs))
            f.write("\n")

def main():
    """
    主函数，读取结构文件并计算原子对距离。
    """
    info = """
    使用说明:
    python get_bond_total.py
    不需要输入任何参数 
    """
    print(info)

    # 从POSCAR文件中读取结构
    struct = Structure.from_file("POSCAR")

    # 计算原子对距离
    pairs_dist = calculate_pair_distances(struct)

    # 打印和保存信息
    print_structure_info(struct)
    print_distance_info(pairs_dist)
    save_distance_info(pairs_dist)

if __name__ == '__main__':
    main()
