#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
脚本名称: cif2vasp.py
用途: 将CIF文件转换为VASP文件格式，并调整z方向的晶格参数。
用法: python cif2vasp_with_z_adjustment.py *.cif
"""

import numpy as np
import sys
import os
from ase.io import read, write
from ase.build import sort

# 获取输入的CIF文件路径
path = sys.argv[1]

# 获取文件名（不带扩展名）
name = os.path.splitext(path)[0]

# 读取CIF文件
structure = read(path)

# 获取晶胞参数
x, y, z = structure.get_cell()

# 获取原子在z方向上的最大最小位置
zmax = structure.get_positions()[:, 2].max()
zmin = structure.get_positions()[:, 2].min()

# 计算新的z方向晶胞参数
z_3 = zmax - zmin
structure.set_cell([x, y, (0.0, 0.0, z_3)])

# 调整原子位置，使其z方向上从0开始
pos = structure.get_positions()
new_pos = []
for i in range(len(pos)):
    x_i = pos[i][0]
    y_i = pos[i][1]
    z_i = pos[i][2] - zmin
    new_pos.append([x_i, y_i, z_i])
new_pos = np.array(new_pos)
structure.set_positions(new_pos)

# 输出VASP格式文件
write('./{}.vasp'.format(name), structure, direct=False, sort=True)
