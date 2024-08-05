# 文件名: get_pairs.py
# 运行方法:
#     python get_pairs.py [COHPstartEnergy] [COHPendEnergy] [element1] [element2] [lower_d_limit] [upper_d_limit]
# 功能描述:
#     指定起始和终止计算的能量，指定两种元素及其间的键长范围，生成Lobster输入文件。
# 注意:
#     Lobster认为你指定的原子对是不具有周期性的。
#     使用pymatgen找到的距离是包含周期性的，把这些原子对输入给Lobster时，它认不出这个距离是周期性的，
#     因此会按照原胞内的距离考虑两个原子的成键。

import os
import sys
from pymatgen.core.structure import Structure

def generate_lobster_input(cohp_start_energy, cohp_end_energy, element1, element2, lower_d, upper_d):
    """
    生成Lobster输入文件。

    参数:
        cohp_start_energy (str): COHP起始能量。
        cohp_end_energy (str): COHP终止能量。
        element1 (str): 元素1。
        element2 (str): 元素2。
        lower_d (float): 原子对之间的最小距离。
        upper_d (float): 原子对之间的最大距离。
    """
    struct = Structure.from_file("POSCAR")

    # 生成lobsterinIsmear_5文件
    with open("lobsterinIsmear_5", "w") as f:
        f.write(f'COHPstartEnergy  {cohp_start_energy}\n')
        f.write(f'COHPendEnergy    {cohp_end_energy}\n')
        f.write('usebasisset pbeVaspFit2015\n')  # 基组（直接使用根据vasp拟合的基组以及默认的基函数）
        for spe in struct.types_of_specie:
            f.write(f'basisfunctions {spe.name}\n')  # 基组（直接使用根据vasp拟合的基组以及默认的基函数）
        f.write(f"cohpGenerator from {lower_d} to {upper_d} type {element1} type {element2}\n")

    # 生成lobsterinIsmear_0文件
    with open("lobsterinIsmear_0", "w") as f:
        f.write(f'COHPstartEnergy  {cohp_start_energy}\n')
        f.write(f'COHPendEnergy    {cohp_end_energy}\n')
        f.write('usebasisset pbeVaspFit2015\n')  # 基组（直接使用根据vasp拟合的基组以及默认的基函数）
        f.write('gaussianSmearingWidth 0.05\n')
        for spe in struct.types_of_specie:
            f.write(f'basisfunctions {spe.name}\n')  # 基组（直接使用根据vasp拟合的基组以及默认的基函数）
        f.write(f"cohpGenerator from {lower_d} to {upper_d} type {element1} type {element2}\n")

def create_output_directory(element1, element2, lower_d, upper_d):
    """
    创建用于存储Lobster输出的目录。

    参数:
        element1 (str): 元素1。
        element2 (str): 元素2。
        lower_d (float): 原子对之间的最小距离。
        upper_d (float): 原子对之间的最大距离。
    """
    dirs = f"{element1}_{element2}_{lower_d}_{upper_d}"
    if not os.path.exists(dirs):
        os.mkdir(dirs)
    return dirs

def print_instructions(dirs):
    """
    打印在运行Lobster前后的指令。

    参数:
        dirs (str): 输出目录名称。
    """
    print("Note: ------------------------------------------------------")
    print("    You need to do it before running Lobster-4.1.0")
    print("    cp lobsterraw lobsterin")
    print("    You need to do it after running Lobster-4.1.0")
    print(f"    cp {{COBICAR.lobster,COHPCAR.lobster,COOPCAR.lobster,ICOBILIST.lobster,ICOHPLIST.lobster,ICOOPLIST.lobster,lobsterin}}  {dirs}")
    print("------------------------------------------------------------")

def main():
    """
    主函数，解析命令行参数并生成Lobster输入文件。
    """
    if len(sys.argv) != 7:
        print("Usage: python get_pairs.py [COHPstartEnergy] [COHPendEnergy] [element1] [element2] [lower_d_limit] [upper_d_limit]")
        sys.exit(1)

    cohp_start_energy = sys.argv[1]
    cohp_end_energy = sys.argv[2]
    element1 = sys.argv[3]
    element2 = sys.argv[4]

    try:
        lower_d = float(sys.argv[5])
        upper_d = float(sys.argv[6])
    except ValueError:
        lower_d = 0.0
        upper_d = 5.0

    generate_lobster_input(cohp_start_energy, cohp_end_energy, element1, element2, lower_d, upper_d)
    dirs = create_output_directory(element1, element2, lower_d, upper_d)
    print_instructions(dirs)

if __name__ == '__main__':
    main()
