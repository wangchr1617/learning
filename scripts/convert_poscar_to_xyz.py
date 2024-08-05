# 文件名: convert_poscar_to_xyz.py
# 运行方法: python convert_poscar_to_xyz.py POSCAR
# 功能描述: 将POSCAR文件转换为XYZ格式的文件

import sys
from ase.io import read, write

def main():
    # 读取命令行参数中的POSCAR文件
    input_file = sys.argv[1]

    # 读取结构信息
    structure = read(input_file)

    # 将结构信息复制为一个超级晶胞
    supercell = structure * (1, 1, 1)  # 这里可以调整超级晶胞的大小

    # 将结构信息写入XYZ文件
    write("model.xyz", supercell)

    # 输出完成提示
    print("转换完成！")

if __name__ == "__main__":
    main()
