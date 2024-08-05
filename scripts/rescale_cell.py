"""
用法: python rescale_cell.py

功能描述:
该脚本读取VASP的POSCAR文件，按指定的缩放因子对晶胞的c轴进行缩放，并生成相应的VASP文件。

参数说明:
- `POSCAR`: 输入的VASP POSCAR文件，必须位于脚本同一目录下。

结果:
- 生成多个VASP文件，文件名格式为`[factor].vasp`，其中`factor`是缩放因子。
"""

from ase.io import read, write

def main():
    # 读取VASP的POSCAR文件
    poscar_file = 'POSCAR'
    original_structure = read(poscar_file)

    # 定义c轴的缩放因子
    scaling_factors = [0.9, 0.92, 0.94, 0.96, 0.98, 1.0]

    # 遍历缩放因子并生成相应的VASP文件
    for factor in scaling_factors:
        modified_structure = original_structure.copy()
        cell = modified_structure.cell.copy()
        
        # 按比例缩放c轴
        cell[2, 2] *= factor
        modified_structure.set_cell(cell, scale_atoms=True)
        
        # 生成输出文件名
        output_filename = f"{factor:.2f}.vasp"
        
        # 写入VASP文件
        write(output_filename, modified_structure, format='vasp', direct=True, vasp5=True)
        print(f"Generated {output_filename}")

    print("所有VASP文件已成功生成。")

if __name__ == "__main__":
    main()
