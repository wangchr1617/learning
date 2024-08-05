# 文件名: filter_xyz.py
# 运行方法: python filter_xyz.py
# 功能描述: 过滤XYZ文件中的帧，删除不符合给定阈值条件的帧，并将结果保存到新的XYZ文件中。

def parse_lattice_line(line):
    """
    解析XYZ文件中晶格参数的行。

    参数:
        line (str): 包含晶格参数的行。

    返回:
        list: 解析后的晶格参数列表。
    """
    parts = line.split('Lattice="')[1].split('"')[0].split()
    lattice_params = [float(x) for x in parts]
    return lattice_params

def parse_force_lines(lines):
    """
    解析XYZ文件中原子力的行。

    参数:
        lines (list): 包含原子力的行列表。

    返回:
        list: 解析后的原子力列表。
    """
    forces = []
    for line in lines:
        force_values = [float(x) for x in line.split()[4:7]]
        forces.append(force_values)
    return forces

def should_delete_frame(lattice_params, forces, force_threshold, lattice_threshold):
    """
    判断是否应该删除某个帧。

    参数:
        lattice_params (list): 帧的晶格参数。
        forces (list): 帧中所有原子的力。
        force_threshold (float): 力的阈值。
        lattice_threshold (float): 晶格参数的阈值。

    返回:
        bool: 是否应该删除该帧。
    """
    # 检查晶格参数是否超过阈值
    if any(abs(param) > lattice_threshold for param in lattice_params):
        return True

    # 检查原子力是否超过阈值
    for force in forces:
        if any(abs(f) > force_threshold for f in force):
            return True

    return False

def process_xyz_file(input_filename, output_filename, force_threshold, lattice_threshold):
    """
    处理XYZ文件，过滤不符合条件的帧。

    参数:
        input_filename (str): 输入的XYZ文件名。
        output_filename (str): 输出的XYZ文件名。
        force_threshold (float): 力的阈值。
        lattice_threshold (float): 晶格参数的阈值。
    """
    with open(input_filename, 'r') as infile, open(output_filename, 'w') as outfile:
        while True:
            line = infile.readline()
            if not line:
                break 
            atom_count = int(line.strip())
            info_line = infile.readline()
            lattice_params = parse_lattice_line(info_line)
            atom_lines = [infile.readline() for _ in range(atom_count)]
            forces = parse_force_lines(atom_lines)
            if not should_delete_frame(lattice_params, forces, force_threshold, lattice_threshold):
                outfile.write(line)
                outfile.write(info_line)
                outfile.writelines(atom_lines)

def main():
    """
    主函数，设置阈值并调用处理函数。
    """
    input_filename = 'dump.xyz'
    output_filename = 'new_dump.xyz'
    force_threshold = 100
    lattice_threshold = 25
    process_xyz_file(input_filename, output_filename, force_threshold, lattice_threshold)
    print(f"Filtered data saved to {output_filename}")

if __name__ == "__main__":
    main()
