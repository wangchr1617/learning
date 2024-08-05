# 文件名: cp2k_pos_simplify.py
# 运行方法: python cp2k_pos_simplify.py <input_file> <interval>
# 功能描述: 从CP2K输出文件中每隔指定间隔提取原子结构数据，并保存到新的XYZ文件中。

import sys
import linecache

def extract_structures(input_file, interval):
    """
    从输入文件中提取结构数据。

    参数:
        input_file (str): 输入文件路径。
        interval (int): 提取间隔。

    返回:
        None
    """
    # 读取输入文件的第一行以获取原子数量
    num_of_atoms = linecache.getline(input_file, 1).strip()
    num_of_atoms = int(num_of_atoms)
    print('Number of Atoms:', num_of_atoms)

    # 计算每个结构在文件中的行数
    lines_per_structure = num_of_atoms + 2

    # 计算文件中的总行数
    with open(input_file, 'r') as file:
        total_lines = sum(1 for _ in file)
    print('Total number of lines:', total_lines)

    # 计算可以提取的结构总数
    total_structures = total_lines // lines_per_structure
    num_of_structures_to_extract = total_structures // interval
    print('Number of structures to extract:', num_of_structures_to_extract)

    # 打开输出文件
    with open('newpos.xyz', 'w') as output:
        # 按指定间隔提取结构
        for i in range(num_of_structures_to_extract):
            for j in range(1, lines_per_structure + 1):
                # 计算行号并提取对应行
                line_number = i * interval * lines_per_structure + j
                line = linecache.getline(input_file, line_number).strip()
                output.write(line + '\n')

def main():
    """
    主函数，解析命令行参数并调用提取函数。
    """
    if len(sys.argv) != 3:
        print("Usage: python repick_xyz.py <input_file> <interval>")
        sys.exit(1)

    input_file = sys.argv[1]
    interval = int(sys.argv[2])

    extract_structures(input_file, interval)

if __name__ == '__main__':
    main()
