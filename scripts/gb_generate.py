# 文件名: gb_generate.py
# 运行方法: python gb_generate.py <o1> <o2> <o3> <sigma> <p1> <p2> <p3>
# 功能描述: 使用aimsgb库生成指定晶界方向和西格玛值的晶界结构，并保存为VAS文件格式。

import sys
from aimsgb import GrainBoundary, Grain

def generate_gb(o_vector, sigma, p_vector, input_file="POSCAR", uc_a=1, uc_b=1):
    """
    生成晶界结构。

    参数:
        o_vector (list): 晶界方向向量。
        sigma (int): 西格玛值。
        p_vector (list): 晶界平移向量。
        input_file (str, optional): 输入结构文件名。默认为"POSCAR"。
        uc_a (int, optional): 晶界结构中晶胞a轴的重复次数。默认为1。
        uc_b (int, optional): 晶界结构中晶胞b轴的重复次数。默认为1。

    返回:
        GrainBoundary: 生成的晶界结构对象。
    """
    # 从输入文件中读取晶粒信息
    s_input = Grain.from_file(input_file)
    
    # 创建晶界结构
    gb = GrainBoundary(o_vector, sigma, p_vector, s_input, uc_a=uc_a, uc_b=uc_b)
    
    return gb

def save_structure_to_file(gb, filename):
    """
    将晶界结构保存为文件。

    参数:
        gb (GrainBoundary): 生成的晶界结构对象。
        filename (str): 输出文件名。
    """
    # 将两个晶粒堆叠为一个晶界结构
    structure = Grain.stack_grains(gb.grain_a, gb.grain_b, direction=gb.direction)
    
    # 保存结构为VAS文件格式
    structure.to(filename, fmt="POSCAR")

def main():
    """
    主函数，解析命令行参数并生成晶界结构。
    """
    if len(sys.argv) != 8:
        print("Usage: python gb_generate.py <o1> <o2> <o3> <sigma> <p1> <p2> <p3>")
        sys.exit(1)

    # 获取命令行参数
    o1, o2, o3 = map(int, sys.argv[1:4])
    sigma = int(sys.argv[4])
    p1, p2, p3 = map(int, sys.argv[5:8])

    # 定义输出文件名
    filename = f"GB_{o1}{o2}{o3}-{sigma}-{p1}{p2}{p3}.vasp"

    # 生成晶界结构
    gb = generate_gb([o1, o2, o3], sigma, [p1, p2, p3])

    # 保存晶界结构到文件
    save_structure_to_file(gb, filename)
    print(f"Grain boundary structure saved to {filename}")

if __name__ == '__main__':
    main()
