#!/usr/bin/env python
# -*- coding: utf-8 -*-

# 文件名: descriptor_analysis_and_visualization.py
# 运行方法: python descriptor_analysis_and_visualization.py
# 功能描述: 读取原子结构数据，计算描述符，使用PCA进行降维，并选择最远点进行可视化。

import numpy as np
import matplotlib.pyplot as plt
from ase.io import read
from pynep.calculate import NEP
from pynep.select import FarthestPointSample
from sklearn.decomposition import PCA

def read_data(filename):
    """
    读取原子结构数据。

    参数:
        filename (str): 原子结构文件路径。

    返回:
        list: 原子结构列表。
    """
    return read(filename, ':')

def calculate_descriptors(atoms, calc):
    """
    计算原子结构的描述符。

    参数:
        atoms (list): 原子结构列表。
        calc (NEP): NEP计算对象。

    返回:
        np.ndarray: 所有原子的描述符矩阵。
        np.ndarray: 描述符的均值。
        np.ndarray: 描述符的标准差。
    """
    descriptors = [calc.get_property('descriptor', atom) for atom in atoms]
    all_descriptors = np.concatenate(descriptors, axis=0)
    descriptor_mean = all_descriptors.mean(axis=0)
    descriptor_std = all_descriptors.std(axis=0)
    return descriptors, descriptor_mean, descriptor_std

def calculate_structure_descriptors(descriptors, descriptor_mean, descriptor_std):
    """
    计算结构的平均描述符，并进行归一化。

    参数:
        descriptors (list): 描述符列表。
        descriptor_mean (np.ndarray): 描述符均值。
        descriptor_std (np.ndarray): 描述符标准差。

    返回:
        np.ndarray: 所有结构的平均描述符矩阵。
        np.ndarray: 所有结构的归一化平均描述符矩阵。
    """
    all_descriptors_mean = np.array(
        [np.mean(descriptor.reshape(1, -1), axis=0) for descriptor in descriptors]
    )
    normalized_descriptors_mean = np.array(
        [np.mean((descriptor - descriptor_mean) / descriptor_std, axis=0) for descriptor in descriptors]
    )
    return all_descriptors_mean, normalized_descriptors_mean

def farthest_point_sample(des, atoms):
    """
    选择最远点样本。

    参数:
        des (np.ndarray): 描述符矩阵。
        atoms (list): 原子结构列表。

    返回:
        list: 选择的原子结构索引。
        list: 未选择的原子结构索引。
    """
    sampler = FarthestPointSample(min_distance=0.05)
    selected_i = sampler.select(des, [], max_select=100)
    unselected_i = [i for i in range(len(atoms)) if i not in selected_i]
    # write('./selected.xyz', [atoms[i] for i in selected_i])
    # write('./unselected.xyz', [atoms[i] for i in unselected_i])
    return selected_i, unselected_i

def perform_pca_and_visualize(all_descriptors_mean, normalized_descriptors_mean):
    """
    使用PCA进行降维并可视化。

    参数:
        all_descriptors_mean (np.ndarray): 所有结构的平均描述符矩阵。
        normalized_descriptors_mean (np.ndarray): 所有结构的归一化平均描述符矩阵。
    """
    # 设置plt的参数
    plt.rcParams['font.size'] = 14
    plt.rcParams['axes.unicode_minus'] = False
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'
    
    fig, axes = plt.subplots(ncols=2, figsize=(7, 3), dpi=140)
    for k, ax in enumerate(axes):
        pca = PCA(n_components=2)
        if k == 0:
            pc = pca.fit_transform(all_descriptors_mean)
            selected_i, _ = farthest_point_sample(all_descriptors_mean, a)
            selected_pc = pca.transform(np.array([all_descriptors_mean[i] for i in selected_i]))
            title = 'Unnormalized'
        else:
            pc = pca.fit_transform(normalized_descriptors_mean)
            selected_i, _ = farthest_point_sample(normalized_descriptors_mean, a)
            selected_pc = pca.transform(np.array([normalized_descriptors_mean[i] for i in selected_i]))
            title = 'Normalized'
        ax.scatter(pc[:, 0], pc[:, 1], s=0.8, c='b', alpha=0.5, label='Raw data')
        ax.scatter(selected_pc[:, 0], selected_pc[:, 1], s=0.8, c='r', alpha=0.5, label='Selected data')
        ax.set_xlabel(f'PCA dimension 0 - Var={pca.explained_variance_ratio_[0]:.2f}')
        ax.set_ylabel(f'PCA dimension 1 - Var={pca.explained_variance_ratio_[1]:.2f}')
        ax.set_title(title)
    ax.legend(frameon=False)
    fig.align_ylabels()
    plt.tight_layout()
    plt.savefig('fps.png', bbox_inches='tight')

if __name__ == '__main__':
    # 读取原子结构数据
    a = read_data('./nepmd.xyz')
    
    # 读取NEP计算对象
    calc = NEP("./nep.txt")
    print(calc)

    # 计算描述符
    descriptors, descriptor_mean, descriptor_std = calculate_descriptors(a, calc)

    # 计算结构描述符
    all_descriptors_mean, normalized_descriptors_mean = calculate_structure_descriptors(descriptors, descriptor_mean, descriptor_std)

    # 执行PCA并可视化
    perform_pca_and_visualize(all_descriptors_mean, normalized_descriptors_mean)
