# 文件名: calculate_distances.py
# 运行方法: python calculate_distances.py <input_file>
# 功能描述: 从XYZ文件中计算每个时间步内金原子(Au)的质心与各个金原子的距离，并将结果保存到文件中。

import sys
import linecache
import numpy as np

def calculate_atom_to_center_distances(file_name, output_file):
    """
    计算金原子(Au)到质心的距离，并将结果写入文件。

    参数:
        file_name (str): 输入文件名。
        output_file (str): 输出文件名。
    """
    # 读取文件的所有行
    with open(file_name, 'r') as f:
        lines = f.readlines()
    print('Total lines in file:', len(lines))

    # 获取原子数
    atom_num = int(linecache.getline(file_name, 1).strip())
    print('Total atom number:', atom_num)

    # 计算迭代次数
    loop_n = len(lines) // (atom_num + 2)
    print('Iteration number:', loop_n)

    # 存储距离的列表
    distances = []

    # 计算每个时间步内的金原子(Au)到质心的距离
    for i in range(loop_n):
        coord1_sum = np.zeros(3)  # 初始化质心坐标累加器
        count_au = 0  # 金原子计数

        # 计算质心
        for j in range(atom_num + 2):
            tmpl = lines[i * (atom_num + 2) + j].split()
            if tmpl[0] == 'Au':
                count_au += 1
                coord1_sum += np.array([float(tmpl[1]), float(tmpl[2]), float(tmpl[3])])

        # 质心坐标
        coord_mass_center1 = coord1_sum / count_au
        c1 = np.array(coord_mass_center1)

        # 计算金原子到质心的距离
        for j in range(atom_num + 2):
            tmpl = lines[i * (atom_num + 2) + j].split()
            if tmpl[0] == 'Au':
                c2 = np.array([float(tmpl[1]), float(tmpl[2]), float(tmpl[3])])
                distances.append(str(np.linalg.norm(c1 - c2)))

    # 将距离写入文件
    with open(output_file, 'w') as f_out:
        f_out.write('\n'.join(distances))

    print(f"The distances between Au atoms and the cluster's mass center have been written to {output_file}")

def main():
    """
    主函数，解析命令行参数并调用计算函数。
    """
    if len(sys.argv) != 2:
        print("Usage: python calculate_distances.py <input_file>")
        sys.exit(1)

    file_name = sys.argv[1]
    output_file = "distance.xyz"

    calculate_atom_to_center_distances(file_name, output_file)

if __name__ == '__main__':
    print('#######  START  #######')
    main()
