#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
用法: python script.py 输入文件.vasp 输出文件.vasp

该脚本读取VASP格式的结构文件，并在特定区域插入Li原子，然后保存修改后的结构到新的文件中。
"""

from ase import Atom
from ase.io import read, write
from ase.build import sort
import sys

def main(filein, fileout):
    # 读取并排序结构
    structure = sort(read(filein))

    # 找到z轴的最大值和最小值
    zmax = structure.get_positions()[:, 2].max()
    zmin = structure.get_positions()[:, 2].min()

    # 找到B原子的位置
    def region(structure, X):
        AtomList = structure.get_chemical_symbols()
        begin = AtomList.index(X)
        for i in range(begin, len(AtomList)):
            if AtomList[i] != X:
                end = i
                break
        else:
            end = len(AtomList)
        return begin, end

    begin, end = region(structure, 'B')

    # 插入Li原子
    for i in range(begin, begin + 4):
        x, y, z = structure.get_positions()[i]
        z1 = zmax + 3
        z2 = zmin - 3
        structure += Atom('Li', (x, y, z1))
        structure += Atom('Li', (x, y, z2))

    # 保存修改后的结构
    write(fileout, structure)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("用法: python script.py 输入文件.vasp 输出文件.vasp")
        sys.exit(1)

    filein = sys.argv[1]
    fileout = sys.argv[2]
    main(filein, fileout)
