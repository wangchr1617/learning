#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
用法: python calculate_bond_angle.py <POSCAR路径> <原子索引1> <原子索引2> <原子索引3>

该脚本加载POSCAR文件并计算指定三个原子的键角。
"""

import sys
from ase.io import read
from ase.geometry import get_angles

def calculate_angle(poscar_path, atom_indices):
    """
    计算给定三个原子之间的键角。

    参数:
        poscar_path: POSCAR文件的路径
        atom_indices: 包含三个原子索引的列表

    返回:
        float，键角（以度为单位）
    """
    # 加载POSCAR文件
    atoms = read(poscar_path)

    # 确保输入的原子索引是有效的
    if all(i < len(atoms) for i in atom_indices):
        # 计算并返回键角
        a, b, c = atom_indices
        v1 = atoms.positions[a] - atoms.positions[b]
        v2 = atoms.positions[c] - atoms.positions[b]
        angle = get_angles(v1, v2)
        return angle[0]
    else:
        raise IndexError("指定的原子索引超出了范围。")

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("使用方法: python calculate_bond_angle.py <POSCAR路径> <原子索引1> <原子索引2> <原子索引3>")
    else:
        poscar_path = sys.argv[1]
        atom_indices = [int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])]
        try:
            angle = calculate_angle(poscar_path, atom_indices)
            print(f"键角: {angle:.2f} 度")
        except IndexError as e:
            print(e)
