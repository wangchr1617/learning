########################################################################################
# Python script to calculate RDF of specified atom and orbital                        ##
# Written by Yafei Jiang                                                              ##
# Email: jiangyafei730@163.com                                                        ##
# Usage: python RDF_analysis.py filename orbital                                      ## 
# Example: python RDF_analysis.py out.Os 5d                                           ##
# Note: the file such as out.Os should be obtained from ADF calculation               ##                                                           ##
########################################################################################

import sys
import subprocess
import numpy as np
import math

def main():
    # 检查命令行参数
    if len(sys.argv) != 3:
        print("Usage: python RDF_analysis.py filename orbital")
        sys.exit(1)

    # 从命令行参数中获取文件名和轨道信息
    file_name = sys.argv[1]
    orbital = sys.argv[2]

    # 读取ADF输出文件
    with open(file_name, 'r') as file_out:
        lines = file_out.readlines()

    # 提取基组信息
    status1, line1 = subprocess.getstatusoutput(f"sed -n '/Valence Basis Sets/=' {file_name}")
    line_basis = int(line1.split('\n')[0])
    Nbasis = int(lines[line_basis - 1].split()[-1])
    
    basis = []  # Zeta values
    mainquantumnumber = []  # Main quantum numbers
    orbitals = []  # Orbital types (S, P, D, F)
    
    for i in range(Nbasis):
        basis.append(float(lines[line_basis + 1 + i].split()[2]))
        mainquantumnumber.append(int(lines[line_basis + 1 + i].split()[0]))
        orbitals.append(lines[line_basis + 1 + i].split()[1])

    # 轨道字典
    orbital_dic = {"S": 0, "P": 1, "D": 2, "F": 3}
    
    # 解析输入的轨道
    orbital_n = int(orbital[0])  # 主量子数
    orbital_ch = orbital[1].upper()  # 角量子数字符
    orbital_l = orbital_dic.get(orbital_ch)  # 角量子数索引

    # 检查角量子数是否有效
    if orbital_l is None:
        print("The orbital you input is not correct.")
        sys.exit(1)

    # 提取系数矩阵
    c, z, N = extract_coefficients(file_name, lines, orbital_ch, orbital_n, orbital_l, basis, mainquantumnumber, orbitals)

    # 计算RDF
    calculate_rdf(file_name, orbital, c, z, N)


def extract_coefficients(file_name, lines, orbital_ch, orbital_n, orbital_l, basis, mainquantumnumber, orbitals):
    """
    提取系数矩阵。

    参数:
    - file_name: 输出文件名。
    - lines: 文件的行内容。
    - orbital_ch: 角量子数字符。
    - orbital_n: 主量子数。
    - orbital_l: 角量子数索引。
    - basis: Zeta值列表。
    - mainquantumnumber: 主量子数列表。
    - orbitals: 轨道类型列表。

    返回:
    - c: 系数矩阵。
    - z: Zeta值。
    - N: 主量子数。
    """
    c = []  # 系数矩阵
    search_pattern = {
        "F": "=== F:xyz ===",
        "D": "=== D:xz ===",
        "P": "=== P:y ===",
        "S": "=== S ==="
    }
    
    status, line = subprocess.getstatusoutput(f"sed -n '/{search_pattern[orbital_ch]}/=' {file_name}")
    
    if status == 0:
        column = orbital_n - orbital_l
        line_coeff = int(line.split('\n')[0])
        index = [i for i in range(len(orbitals)) if orbitals[i] == orbital_ch]
        z = [basis[i] for i in index]
        N = [mainquantumnumber[i] for i in index]

        # 读取系数
        for i in range(len(z)):
            c_value = float(lines[line_coeff + 6 + i].split()[column])
            c.append(c_value)
    else:
        print(f"Failed to find coefficients for orbital {orbital_ch}.")
        sys.exit(1)

    return c, z, N


def calculate_rdf(file_name, orbital, c, z, N):
    """
    计算径向分布函数 (RDF)。

    参数:
    - file_name: 输入文件名。
    - orbital: 输入的轨道。
    - c: 系数矩阵。
    - z: Zeta值。
    - N: 主量子数。
    """
    M = len(z)
    
    # 原子轨道 (AO) 归一化系数
    b = np.zeros(M)
    for i in range(M):
        b[i] = (2 * z[i]) ** (N[i] + 0.5) / (math.factorial(2 * N[i])) ** 0.5
    
    # 分子轨道 (MO) 归一化系数
    V = 0
    for i in range(M):
        for j in range(M):
            V += math.factorial(N[i] + N[j]) * b[i] * b[j] * c[i] * c[j] / (z[i] + z[j]) ** (N[i] + N[j] + 1)

    # 计算RDF
    grid = 800
    D_r = np.zeros(grid)
    r = np.zeros(grid)
    
    for m in range(grid):
        r[m] = m / 100  # 径向变量
        sumsum = 0
        for p in range(M):
            for q in range(M):
                sumpq = (1 / V) * r[m] ** 2 * c[p] * c[q] * b[p] * b[q] * (r[m] ** (N[p] - 1) * np.exp(-z[p] * r[m])) * (r[m] ** (N[q] - 1) * np.exp(-z[q] * r[m]))
                sumsum += sumpq
        D_r[m] = sumsum

    # 将结果写入文件
    bohr2A = 0.529
    output_filename = f"{file_name.split('.')[1]}_{orbital}-RDF.dat"
    
    with open(output_filename, 'w') as f:
        for i in range(len(r)):
            f.write(f"{r[i] * bohr2A:^10.6f}{D_r[i]:^10.6f}\n")

    print("RDF calculation is finished.")
    print(f"Results saved to {output_filename}")


if __name__ == "__main__":
    main()
