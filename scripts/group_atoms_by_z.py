# 文件名: group_atoms_by_z.py
# 运行方法: python group_atoms_by_z.py
# 功能描述: 根据原子在z轴上的坐标将其分组，并将带有组属性的原子信息写入新的XYZ文件。

def group_atoms(atoms, ranges):
    """
    根据z坐标将原子分组。

    参数:
        atoms (list): 原子列表，每个原子是一个字典，包含元素类型和坐标。
        ranges (list): 分组范围列表，每个范围是一个三元组(y_min, y_max, group_property)。
    """
    for atom in atoms:
        # 取得原子的z坐标
        z = atom['z']
        # 遍历定义的分组范围
        for z_min, z_max, group_property in ranges:
            # 检查z坐标是否在当前范围内
            if z_min <= z < z_max:
                # 如果在范围内，分配对应的组属性
                atom['property'] = group_property
                break

def read_xyz_file(filename):
    """
    读取XYZ文件，并返回原子列表和晶格属性行。

    参数:
        filename (str): 输入XYZ文件的路径。

    返回:
        tuple: 包含原子列表和晶格属性行的元组。
    """
    with open(filename, 'r') as f:
        lines = f.readlines()

    # 读取第二行的晶格属性信息
    lattice_properties_line = lines[1].strip()
    # 添加自定义组属性说明
    lattice_properties_line += ":group:I:1\n"

    # 初始化原子列表
    atoms = []
    # 从第三行开始解析原子数据
    for line in lines[2:]:
        atom_data = line.split()
        atom = {
            'element': atom_data[0],
            'x': float(atom_data[1]),
            'y': float(atom_data[2]),
            'z': float(atom_data[3]),
            'property': None
        }
        atoms.append(atom)

    return atoms, lattice_properties_line

def write_xyz_file(filename, atoms, lattice_properties_line):
    """
    将带有组属性的原子信息写入新的XYZ文件。

    参数:
        filename (str): 输出XYZ文件的路径。
        atoms (list): 带有组属性的原子列表。
        lattice_properties_line (str): 晶格属性行。
    """
    with open(filename, 'w') as f:
        # 写入原子数量
        f.write(f"{len(atoms)}\n")
        # 写入晶格属性行
        f.write(lattice_properties_line)
        # 写入每个原子的详细信息
        for atom in atoms:
            f.write(f"{atom['element']} {atom['x']} {atom['y']} {atom['z']} {atom['property']}\n")

def main():
    """
    主函数，读取原子文件，分组，并输出带有组属性的新原子文件。
    """
    # 定义分组范围，格式为(z_min, z_max, group_property)
    ranges = [(-1, 64, 0), (64, 127, 1)]

    # 读取原子数据和晶格属性行
    atoms, lattice_properties_line = read_xyz_file("model.xyz")

    # 对原子进行分组
    group_atoms(atoms, ranges)

    # 写入新的XYZ文件
    write_xyz_file("model_new.xyz", atoms, lattice_properties_line)

    print("New XYZ file has been generated.")

if __name__ == "__main__":
    main()
