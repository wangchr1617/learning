"""
文件名: create_chemiscope_input.py
运行方法: python create_chemiscope_input.py -a <cutoff_value> <ase_dataset.db> <descriptors.npy> <energies.npy>
功能描述: 该脚本从 ASE 数据库中读取结构数据，并使用描述符进行 PCA 降维，生成适用于 Chemiscope 的输入文件。
"""

import argparse
import numpy as np
import chemiscope as cs
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from ase.io import read

def main():
    # 创建解析器对象
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)

    # 添加命令行参数
    parser.add_argument('dataset', help='ASE 数据集的路径')
    parser.add_argument('descriptors', help='数据集对应的描述符文件路径')
    parser.add_argument('energies', help='数据集对应的能量文件路径')
    parser.add_argument('-a', '--atomic', type=float, help='指定局部环境的截断距离')

    # 解析命令行参数
    args = parser.parse_args()

    print('读取数据中...')
    # 读取 ASE 数据集
    frames = read(args.dataset, ':')
    # 读取描述符和能量数据
    descriptors = np.load(args.descriptors)
    energies = np.load(args.energies)

    print('执行 PCA...')
    # 标准化描述符数据
    scaler = StandardScaler()
    descriptors_scaled = scaler.fit_transform(descriptors)

    # PCA 降维到三维
    pca = PCA(n_components=3, svd_solver='full')
    descriptors_pc = pca.fit_transform(descriptors_scaled)

    print('创建 Chemiscope 输入文件...')
    # 判断目标是结构还是原子
    target = 'structure' if args.atomic is None else 'atom'

    # 创建属性字典
    properties = {
        'PCA': {
            'target': target,
            'values': descriptors_pc,
            'description': '描述符的 PCA 降维结果',
        },
        'energies': {
            'target': target,
            'values': energies,
            'units': 'eV',
            'description': '结构的能量（每个原子）',
        },
    }

    # 生成 Chemiscope 输入文件
    # Chemiscope: https://chemiscope.org/
    if target == 'structure':
        cs.write_input(
            path='chemiscope.json.gz',
            frames=frames,
            properties=properties
        )
    else:
        cs.write_input(
            path='chemiscope.json.gz',
            frames=frames,
            properties=properties,
            environments=cs.all_atomic_environments(frames, args.atomic)
        )

    print('完成！请将 chemiscope.json.gz 上传到 https://chemiscope.org/ 查看。')

if __name__ == '__main__':
    main()
