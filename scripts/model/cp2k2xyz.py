"""
文件名: cp2k2xyz.py
运行方法: python cp2k2xyz.py <目录路径>
功能描述: 该脚本从指定目录读取 CP2K 输出文件，并将数据转换为 XYZ 格式。
"""

import sys
import os

def get_cp2k_filenames(directory):
    """
    从指定目录获取 CP2K 输出文件名，包括位置、力和晶胞文件。

    参数:
    directory (str): 包含 CP2K 输出文件的目录路径。

    返回:
    tuple: 包含位置文件、力文件和晶胞文件名的元组。

    抛出:
    FileNotFoundError: 如果没有找到所有所需的文件。
    """
    position_file = None
    force_file = None
    cell_file = None

    # 遍历目录中的文件
    for filename in os.listdir(directory):
        if filename.endswith(".xyz"):
            if "-pos-" in filename:
                position_file = filename
            elif "-frc-" in filename:
                force_file = filename

        elif filename.endswith(".cell"):
            cell_file = filename

    # 检查是否找到所有必需的文件
    if position_file is None or force_file is None or cell_file is None:
        print(f"请检查目录 {directory} 中的文件。")
        raise FileNotFoundError("找不到必要的文件: 位置文件, 力文件或晶胞文件缺失。")

    return position_file, force_file, cell_file

def convert_cp2k_to_xyz(directory, filenames=None, output_filename=None):
    """
    将 CP2K 输出文件转换为 XYZ 格式。

    参数:
    directory (str): 包含 CP2K 输出文件的目录路径。
    filenames (tuple): 包含位置文件、力文件和晶胞文件名的元组。
    output_filename (str): 输出 XYZ 文件的名称。
    """
    if filenames is None:
        # 获取默认文件名
        position_file, force_file, cell_file = get_cp2k_filenames(directory)
    else:
        position_file, force_file, cell_file = filenames

    # 创建输出目录
    os.makedirs("CP2K2XYZ", exist_ok=True)
    if output_filename is None:
        output_filename = "merged_cp2k.xyz"

    # 组合文件路径
    position_path = os.path.join(directory, position_file)
    force_path = os.path.join(directory, force_file)
    cell_path = os.path.join(directory, cell_file)
    output_path = os.path.join("CP2K2XYZ", output_filename)
    
    print("输入文件:", position_path, force_path, cell_path)
    print("输出文件:", output_path)

    # 打开输入和输出文件
    with open(position_path, 'r') as pf, open(force_path, 'r') as ff, \
         open(cell_path, 'r') as cf, open(output_path, 'w') as of:
        cf.readline()  # 跳过晶胞文件的头部行

        while True:
            # 读取头部信息
            position_header = pf.readline()
            force_header = ff.readline()

            if not position_header or not force_header:
                break  # 如果到达文件结尾，则退出循环

            # 解析原子数量
            num_atoms = int(position_header.strip().split()[0])
            of.write(f"{num_atoms}\n")
            force_info_line = ff.readline()
            pf.readline()  # 跳过位置文件的对应行

            # 计算能量并转换单位
            energy = float(force_info_line.strip().split("E =")[-1]) * 27.211386245988

            # 读取并处理晶胞信息
            cell_line = cf.readline().strip().split()
            lattice = " ".join(cell_line[2:11])  # 读取 Ax, Ay, Az, Bx, By, Bz, Cx, Cy, Cz 列

            # 写入输出文件的头部信息
            of.write(f"energy={energy:.10f} config_type=cp2k2xyz pbc=\"T T T\" ")
            of.write(f"Lattice=\"{lattice}\" Properties=species:S:1:pos:R:3:force:R:3\n")

            # 读取并写入每个原子的位置信息和力信息
            for _ in range(num_atoms):
                position_line = pf.readline().strip().split()
                force_line = ff.readline().strip().split()

                if len(position_line) < 4 or len(force_line) < 4:
                    break  # 如果数据格式不正确，则退出

                # 将力转换为所需的单位
                force_x = float(force_line[1]) * 51.42206747632590000
                force_y = float(force_line[2]) * 51.42206747632590000
                force_z = float(force_line[3]) * 51.42206747632590000

                # 写入原子信息
                of.write(f"{position_line[0]} {position_line[1]} {position_line[2]} {position_line[3]} ")
                of.write(f"{force_x:.10f} {force_y:.10f} {force_z:.10f}\n")

if __name__ == "__main__":
    # 检查是否提供了目录路径参数
    if len(sys.argv) != 2:
        print("用法: python cp2k_to_xyz_converter.py <目录路径>")
        sys.exit(1)

    # 获取目录路径
    directory_path = sys.argv[1]
    convert_cp2k_to_xyz(directory_path)
