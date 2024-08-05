# 文件名: lattice_transformation_calculator.py
# 运行方法: python lattice_transformation_calculator.py
# 功能描述: 读取POSCAR文件中的晶格参数，并计算原始单元格和基元单元格之间的转换矩阵。

import numpy as np

def read_poscar_lattice(filename):
    """
    读取POSCAR文件中的晶格参数。

    参数:
        filename (str): POSCAR文件的路径。

    返回:
        np.ndarray: 3x3的晶格矩阵。
    """
    lattice = np.zeros((3, 3))
    try:
        with open(filename, 'r') as file:
            file.readline()  # 跳过第一行注释
            scale_factor = float(file.readline().strip())  # 读取缩放因子
            # 读取晶格参数
            for i in range(3):
                line = file.readline().strip().split()
                lattice[i] = [float(x) * scale_factor for x in line[:3]]
    except IOError:
        print(f"Error: Could not read file {filename}")
    except ValueError:
        print("Error: File format is not correct")
    return lattice

def main():
    """
    主函数，读取晶格参数并计算转换矩阵。
    """
    # 读取原始单元格晶格参数
    unit = np.array(read_poscar_lattice("./POSCAR"))
    print("原始单元格晶格矩阵：")
    print(unit)

    # 读取基元单元格晶格参数
    prim = np.array(read_poscar_lattice("./PRIMCELL.vasp"))
    print("基元单元格晶格矩阵：")
    print(prim)

    # 计算原始单元格晶格矩阵的逆矩阵
    unit_inv = np.linalg.inv(unit)

    # 计算转换矩阵
    result = np.dot(unit_inv, prim)
    print("转换矩阵：")
    print(result)

if __name__ == '__main__':
    main()
