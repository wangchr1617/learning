# 文件名: pos2xyz.py
# 运行方法: python pos2xyz.py POSCAR [--sort_axis x|y|z] [--order asc|desc] [--sx Sx] [--sy Sy] [--sz Sz]
# 功能描述: 将POSCAR文件转换为XYZ格式的文件，并可以选择按指定轴排序输出，控制超胞的复制倍数和排序顺序

import sys
import argparse
from ase.io import read, write
from ase.build import sort

def main():
    # 创建解析器
    parser = argparse.ArgumentParser(description="将POSCAR文件转换为XYZ格式的文件，并可以选择按指定轴排序输出，控制超胞的复制倍数和排序顺序")
    parser.add_argument("input_file", help="输入的POSCAR文件")
    parser.add_argument("--sort_axis", choices=['x', 'y', 'z'], help="按指定轴排序输出XYZ文件")
    parser.add_argument("--order", choices=['asc', 'desc'], default='asc', help="排序顺序，asc为升序，desc为降序，默认为asc")
    parser.add_argument("--sx", type=int, default=1, help="超胞在x方向的复制倍数，默认为1")
    parser.add_argument("--sy", type=int, default=1, help="超胞在y方向的复制倍数，默认为1")
    parser.add_argument("--sz", type=int, default=1, help="超胞在z方向的复制倍数，默认为1")

    # 解析参数
    args = parser.parse_args()

    # 读取结构信息
    structure = read(args.input_file)

    # 将结构信息复制为一个超胞
    supercell = structure * (args.sx, args.sy, args.sz)

    # 如果需要按指定轴排序
    if args.sort_axis:
        sorted_supercell = sort(supercell, tags=supercell.positions[:, 'xyz'.index(args.sort_axis)])
        if args.order == 'desc':
            sorted_supercell = sort(supercell, tags=-supercell.positions[:, 'xyz'.index(args.sort_axis)])
        supercell = sorted_supercell

    # 将结构信息写入XYZ文件
    write("model.xyz", supercell)

    # 输出完成提示
    print("转换完成！")

if __name__ == "__main__":
    main()
