"""
用法: python combine_vasp_structures.py A.vasp B.vasp

功能描述:
该脚本用于合并两个VASP文件中的结构，并将它们保存为一个新的VASP文件。
结构A和结构B在z方向上会以指定的间隔连接起来。

参数说明:
- A.vasp: 第一个VASP文件的路径。
- B.vasp: 第二个VASP文件的路径。

结果:
生成的合并结构将保存为新的VASP文件，文件名为AB.vasp。
"""

from ase.io import read, write
from ase.build import sort
import sys
import numpy as np

def main():
    # 获取命令行参数
    file_A = sys.argv[1]  # './A.vasp'
    file_B = sys.argv[2]  # './B.vasp'
    
    # 读取并排序结构
    structure_A = sort(read(file_A))
    structure_B = sort(read(file_B))
    
    # 获取原子晶胞
    cell_A = structure_A.cell
    cell_B = structure_B.cell

    # 设置合并结构的间隔
    spacing = 2.0
    
    # 更新结构A和结构B的晶胞以匹配xy方向的最大值
    new_cell_A = np.array(cell_A)
    new_cell_A[0:2] = np.maximum(cell_A[0:2], cell_B[0:2])
    new_cell_B = np.array(cell_B)
    new_cell_B[0:2] = np.maximum(cell_A[0:2], cell_B[0:2])

    # 创建新的合并晶胞
    new_cell_C = np.zeros((3, 3))
    new_cell_C[0:2] = np.maximum(cell_A[0:2], cell_B[0:2])
    new_cell_C[2] = cell_A[2] + cell_B[2] + [0, 0, spacing]

    # 更新结构的晶胞，不缩放原子
    structure_A.set_cell(new_cell_A, scale_atoms=False)
    structure_B.set_cell(new_cell_B, scale_atoms=False)

    # 计算结构B的z方向平移量
    zmax_A = structure_A.get_positions()[:, 2].max()
    zmin_B = structure_B.get_positions()[:, 2].min()
    z_translation = zmax_A - zmin_B + spacing

    # 平移结构B
    structure_B.translate((0, 0, z_translation))

    # 合并结构A和结构B
    structure_C = structure_A + structure_B
    structure_C.set_cell(new_cell_C, scale_atoms=False)

    # 保存合并后的结构为新的VASP文件
    write('AB.vasp', structure_C)
    print("合并结构已保存为: AB.vasp")

if __name__ == "__main__":
    main()
