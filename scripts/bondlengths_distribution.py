#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
用法: python bondlengths_distribution.py <结构文件路径> <输出文件标识>

该脚本读取结构文件（XYZ格式），计算所有键长，并生成键长分布的直方图。
"""

import sys
import matplotlib.pyplot as plt
from ase.data import covalent_radii
from ase.io import read
from ase.neighborlist import NeighborList

# 设置画图格式
plt.rcParams['font.size'] = 14
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'

def get_bondpairs(atoms, radius=1.1):
    """
    获取原子对及其距离

    参数:
        atoms: ASE原子对象
        radius: 用于计算邻居列表的半径比例

    返回:
        bondpairs: 包含原子对及其距离的列表
    """
    cutoffs = radius * covalent_radii[atoms.numbers]
    nl = NeighborList(cutoffs=cutoffs, self_interaction=False)
    nl.update(atoms)
    bondpairs = []
    for a in range(len(atoms)):
        indices, _ = nl.get_neighbors(a)
        for b in indices:
            dis = atoms.get_distance(a, b, mic=True)
            bondpairs.append((a, b, dis))
    return bondpairs

def main(structure_file, idx, radius=1.1):
    """
    主函数，读取结构文件并生成键长分布直方图

    参数:
        structure_file: 结构文件路径
        idx: 输出文件标识
        radius: 用于计算邻居列表的半径比例
    """
    traj = read(structure_file, index=':')
    bondlengths = []
    for atoms in traj:
        bondpairs = get_bondpairs(atoms, radius)
        bondlengths += [pair[2] for pair in bondpairs]

    fig, ax = plt.subplots()
    ax.hist(bondlengths, bins=50, density=True, alpha=0.5, color='r', edgecolor='black')
    ax.set_xlabel('Bond Length')
    ax.set_ylabel('Frequency')
    plt.savefig(f'./bondlengths_{idx}.png')

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("用法: python bondlengths_distribution.py <结构文件路径> <输出文件标识>")
    else:
        structure_file = sys.argv[1]
        idx = sys.argv[2]
        main(structure_file, idx)
