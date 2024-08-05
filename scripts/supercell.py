"""
用法: python supercell.py [input_file] [sx] [sy] [sz]

功能描述:
该脚本读取指定的VASP文件，并生成超胞结构。超胞结构的大小由输入参数[sx, sy, sz]指定。

参数说明:
- [input_file]: VASP格式的输入文件路径。
- [sx]: x方向的超胞倍数。
- [sy]: y方向的超胞倍数。
- [sz]: z方向的超胞倍数。

结果:
生成的超胞结构将保存为新的VASP文件，文件名格式为 [input_file]_SC[sx][sy][sz].vasp。
"""

import sys
import os
from ase.io import read, write
from ase.build import sort

def main():
    # 获取命令行参数
    path = sys.argv[1]
    sx = int(sys.argv[2])
    sy = int(sys.argv[3])
    sz = int(sys.argv[4])

    # 从文件路径中提取文件名（不包括扩展名）
    name = os.path.splitext(path)[0]

    # 读取结构并生成超胞
    structure = read(path) * (sx, sy, sz)

    # 保存超胞结构为新的VASP文件，排序原子并使用direct坐标
    output_filename = f'./{name}_SC{sx}{sy}{sz}.vasp'
    write(output_filename, sort(structure), direct=True)
    print(f"超胞结构已保存为: {output_filename}")

if __name__ == "__main__":
    main()
