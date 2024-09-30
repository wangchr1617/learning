"""
脚本名称: cif2vasp.py
用途: 将CIF文件转换为VASP文件格式。
用法: python cif2vasp.py mp-2612.cif
"""

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

# 将结构按原子类型排序并写入VASP格式文件
write('./{}.vasp'.format(name), sort(structure), direct=True)
