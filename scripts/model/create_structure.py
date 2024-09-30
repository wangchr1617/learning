# 文件名: create_structure.py
# 运行方法: python create_structure.py
# 功能描述: 使用ASE库创建一个晶体结构并将其保存为VASP文件格式。

from ase import Atoms
from ase.build import sort 
from ase.io import write

def create_gete_structure():
    """
    创建GeTe晶体结构并保存为VASP文件。
    """
    # 定义晶体结构参数
    lattice_para = 6.015  # 晶格参数
    angle = 89.954  # 晶体角度

    # 定义原子结构
    atoms = Atoms('Ge4Te4', 
                  cell=[lattice_para, lattice_para, lattice_para, angle, angle, angle], 
                  pbc=(1, 1, 1),  # 周期性边界条件
                  scaled_positions=[(0.0, 0.0, 0.0), (0.0, 0.5, 0.5), (0.5, 0.0, 0.5), (0.5, 0.5, 0.0),
                                    (0.5, 0.0, 0.0), (0.5, 0.5, 0.5), (0.0, 0.0, 0.5), (0.0, 0.5, 0.0)])
    
    # 将结构排序并写入VAS文件
    write("GeTe_cubic.vasp", sort(atoms), direct=True)

def main():
    """
    主函数，调用创建结构函数。
    """
    create_gete_structure()

if __name__ == '__main__':
    main()
