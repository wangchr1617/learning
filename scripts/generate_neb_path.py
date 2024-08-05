# 文件名: generate_neb_path.py
# 运行方法: python generate_neb_path.py <initial_POSCAR> <final_POSCAR> <n_images>
# 功能描述: 使用IDPP方法生成用于NEB计算的插值路径。

import os
import sys
from pymatgen.core import Structure
from pymatgen.analysis.diffusion.neb.pathfinder import IDPPSolver

def generate_neb_path(init_file, final_file, n_images):
    """
    使用IDPP方法生成插值路径，并将每个图像保存为POSCAR文件。

    参数:
        init_file (str): 初始结构文件（POSCAR）的路径。
        final_file (str): 终止结构文件（POSCAR）的路径。
        n_images (int): 图像数量。
    """
    # 读取初始和终止结构
    init_struct = Structure.from_file(init_file, False)
    final_struct = Structure.from_file(final_file, False)

    # 使用IDPPSolver生成插值路径
    idpp_solver = IDPPSolver.from_endpoints(endpoints=[init_struct, final_struct], nimages=n_images, sort_tol=1.0)
    new_path = idpp_solver.run(maxiter=5000, tol=1e-5, gtol=1e-3, step_size=0.05, max_disp=0.05, spring_const=5.0)

    # 保存每个插值结构为POSCAR文件
    for i, structure in enumerate(new_path):
        image_dir = f'{i:02d}'
        if not os.path.exists(image_dir):
            os.makedirs(image_dir)
        poscar_file = os.path.join(image_dir, 'POSCAR')
        structure.to(fmt="poscar", filename=poscar_file)

    print("Improved interpolation of NEB initial guess has been generated. BYE.")

def main():
    """
    主函数，解析命令行参数并调用生成插值路径的函数。
    """
    # 重定向标准输出
    sys.stdout = open(os.devnull, 'w')

    if len(sys.argv) < 4:
        raise SystemError('Syntax Error! Run as: python generate_neb_path.py <initial_POSCAR> <final_POSCAR> <n_images>')

    # 获取命令行参数
    init_file = sys.argv[1]
    final_file = sys.argv[2]
    n_images = int(sys.argv[3])

    # 生成NEB路径
    generate_neb_path(init_file, final_file, n_images)

    # 恢复标准输出
    sys.stdout = sys.__stdout__

if __name__ == '__main__':
    main()
