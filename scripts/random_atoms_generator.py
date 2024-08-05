#!/usr/bin/env python

"""
Young-Jae Choi, POSTECH, Korea, Rep. of.
Inquiries: ssrokyz@gmail.com

用法: python random_atoms_generator.py

功能描述:
该脚本基于给定的骨架结构生成具有随机位置的原子图像。可以指定原子的数量、固定的原子位置、
真空层、应变等多个参数。

依赖库:
- numpy
- ase
"""

import numpy as np
from copy import deepcopy
from ase.atoms import symbols2numbers as s2n
from ase.data import covalent_radii
from ase.build import make_supercell
from ase.optimize.precon.neighbors import estimate_nearest_neighbour_distance as rNN
from time import time

def count2list(dict_in):
    """
    将元素计数的字典转换为列表。

    参数:
    - dict_in: 元素计数字典。

    返回:
    - list_out: 包含重复元素的列表。
    """
    list_out = []
    for key, itera in dict_in.items():
        list_out.extend([key] * itera)
    return list_out

def list2count(list_inp):
    """
    将元素列表转换为计数字典。

    参数:
    - list_inp: 输入的元素列表。

    返回:
    - dict_out: 元素计数字典。
    """
    dict_out = {}
    for element in list_inp:
        if element in dict_out:
            dict_out[element] += 1
        else:
            dict_out[element] = 1
    return dict_out

def covalent_expect(input):
    """
    计算系统的共价键期望值。

    参数:
    - input: 输入可以是字典或列表。{'Si':1, 'H':2} 或 ['Si', 'H', 'H']

    返回:
    - expect_value: 共价键期望值。
    """
    if isinstance(input, list):
        num_spec_dict = list2count(input)
    elif isinstance(input, dict):
        num_spec_dict = input
    else:
        raise TypeError("input is not list nor dict")
    
    tot_num = np.sum(list(num_spec_dict.values()))
    r_sum = sum(covalent_radii[s2n([key])[0]] * value for key, value in num_spec_dict.items())
    
    return r_sum / tot_num

def random_atoms_gen(
    backbone,
    num_spec_dict = None,
    fix_ind_dict  = None,
    pin_the_fixed = False,
    cutoff_radi   = None,
    cutoff_frac   = None,
    random_radi   = None,
    random_frac   = None,
    strain        = None,
    strain_ratio  = [1., 1., 1.],
    vacuum        = None,
    vacuum_ratio  = None,
    sort_atoms    = False,
    max_trial_sec = 5,
    log           = True,
    ):
    """ 
    RAG生成基于骨架结构的随机位置图像。

    参数说明:
    - backbone: ASE原子对象，骨架结构。
    - num_spec_dict: 每种原子的数量字典。如果为None，则与骨架结构相同。"V"表示空位。
    - fix_ind_dict: 每种原子固定的索引字典或列表。如果所有原子位置都必须打乱，则设置为None。
    - pin_the_fixed: 如果为True，则固定原子不偏离其晶格点。
    - cutoff_radi: 两个原子之间最小距离的截断半径。
    - cutoff_frac: 截断半径的比例值，以期望共价半径为基础。
    - random_radi: 最大随机偏移矢量。
    - random_frac: 随机半径的比例值，以RDF的第一个峰值为基础。
    - strain: 应变值列表。
    - strain_ratio: 应变比例列表。
    - vacuum: 真空层调整值列表。
    - vacuum_ratio: 真空层调整比例列表。
    - sort_atoms: 是否按化学序号对原子排序。
    - max_trial_sec: 生成新原子位置的最大尝试时间。
    - log: 是否输出日志信息。
    """
    
    backbone = backbone.copy()
    
    # 获取原子种类列表和数量字典
    if num_spec_dict is None:
        spec_list = backbone.get_chemical_symbols()
        num_spec_dict = list2count(spec_list)
        num_vacancy = 0
    else:
        num_vacancy = num_spec_dict.pop('V', 0)
        spec_list = count2list(num_spec_dict)

    # 获取固定原子信息
    num_fix_dict = {}
    num_shffl_spec_dict = deepcopy(num_spec_dict)
    if isinstance(fix_ind_dict, (list, np.ndarray)):
        fix_ind_arr = np.array(deepcopy(fix_ind_dict))
        fixed_atoms = backbone.copy()[fix_ind_arr]
        fixed_spec = np.array(fixed_atoms.get_chemical_symbols())
        fix_ind_dict = {spec: list(fix_ind_arr[fixed_spec == spec]) for spec in np.unique(fixed_spec)}
    if isinstance(fix_ind_dict, dict): 
        fix_ind_dict = deepcopy(fix_ind_dict)
        for key, value in fix_ind_dict.items():
            num_fix_dict[key] = len(value)
            num_shffl_spec_dict[key] -= len(value)
            fix_ind_dict[key] = np.array(value, dtype=np.int32).tolist()
            if key == 'V' and num_vacancy == 0:
                raise ValueError('fix_ind_dict cannot have "V" if num_spec_dict does not have "V".')
        fix_ind_dict.setdefault('V', [])
        num_fix_dict.setdefault('V', 0)
        num_shffl_spec_dict.setdefault('V', 0)
    elif fix_ind_dict is None:
        fix_ind_dict = {key: [] for key in num_spec_dict}
        fix_ind_dict['V'] = []
        num_fix_dict = {key: 0 for key in num_spec_dict}
        num_fix_dict['V'] = 0
    else:
        raise ValueError('Unknown type of fix_ind_dict.')

    # 获取需要打乱的原子列表
    shffl_spec_list = count2list(num_shffl_spec_dict)

    # 计算共价键长度期望值
    coval_expect = covalent_expect(spec_list)
                
    ## 调整晶胞应变
    if strain_ratio is not None and strain is not None:
        raise ValueError("strain_ratio & strain parameters provided simultaneously. Just provide one.")
    if strain is not None:
        strain = np.array(strain)
        if strain.shape != (3,):
            raise ValueError("Something is wrong with strain parameter. Please check.")
        norm = np.linalg.norm(backbone.cell, axis=1)
        strain_ratio = strain / norm + 1
    if strain_ratio is not None:
        strain_ratio = np.array(strain_ratio)
        if strain_ratio.shape != (3,):
            raise ValueError("Something is wrong with strain_ratio parameter. Please check.")
        backbone.set_cell(
            backbone.cell * np.expand_dims(strain_ratio, axis=1),
            scale_atoms=True,
        )
    if strain_ratio is None and strain is None:
        strain_ratio = [1., 1., 1.]
        backbone.set_cell(
            backbone.cell * np.expand_dims(strain_ratio, axis=1),
            scale_atoms=True,
        )

    ## 调整真空层
    if vacuum_ratio is not None and vacuum is not None:
        raise ValueError("vacuum_ratio & vacuum parameters provided simultaneously. Just provide one.")
    if vacuum is not None:
        vacuum = np.array(vacuum)
        if vacuum.shape != (3,):
            raise ValueError("Something is wrong with vacuum parameter. Please check.")
        norm = np.linalg.norm(backbone.cell, axis=1)
        vacuum_ratio = vacuum / norm + 1
    if vacuum_ratio is not None:
        vacuum_ratio = np.array(vacuum_ratio)
        if vacuum_ratio.shape != (3,):
            raise ValueError("Something is wrong with vacuum_ratio parameter. Please check.")
        backbone.set_cell(
            backbone.cell * np.expand_dims(vacuum_ratio, axis=1),
            scale_atoms=False,
        )
    if vacuum_ratio is None and vacuum is None:
        vacuum_ratio = [1., 1., 1.]
        backbone.set_cell(
            backbone.cell * np.expand_dims(vacuum_ratio, axis=1),
            scale_atoms=True,
        )

    ## 确定截断半径
    if cutoff_radi is not None and cutoff_frac is not None:
        raise ValueError("cutoff_radi & cutoff_frac parameters provided simultaneously. Just provide one.")
    if cutoff_radi is not None:
        cutoff_r = cutoff_radi
    elif cutoff_frac is not None:
        cutoff_r = coval_expect * 2 * cutoff_frac
    else:
        cutoff_r = 0.

    ## 获取随机调整半径
    if random_frac is not None and random_radi is None:
        supercell = make_supercell(backbone, [[2, 0, 0], [0, 2, 0], [0, 0, 2]])
        rdf_1st_peak = rNN(supercell)
        ran_radi = rdf_1st_peak / 2 * random_frac
        if log:
            print("")
            print("********* Please check carefully !!!! ***********".center(120))
            print(f"RDF 1st peak / 2 == {rdf_1st_peak/2:.2f}".center(120))
            print(f"Positional deviation degree == {random_frac:.2f}".center(120))
            print(f"==> Random deviation radius == {ran_radi:.2f}".center(120))
            print(f"it is {ran_radi / coval_expect * 100:.2f} % of covalent-bond-length expectation value.".center(120))
            print("")
            print(f"C.f. ) covalent-bond-length expectation value == {coval_expect:.2f}".center(120))
            print(f"C.f. ) cutoff radius == {cutoff_r:.2f}".center(120))
            print(f"C.f. ) cutoff radius / covalent bond expectation *2 == {cutoff_r / coval_expect / 2 * 100:.2f} %".center(120))
            print("")
    elif random_radi is not None and random_frac is None:
        ran_radi = float(random_radi)
    else:
        raise ValueError('Check random_radi or random_frac parameters.')

    ## 主程序
    if num_vacancy != 0:
        # 选择空位索引
        vacancy_ind = np.random.permutation(
            np.setdiff1d(
                range(len(backbone)),
                np.concatenate(list(fix_ind_dict.values())),
                True,
            ),
        )[:num_vacancy - num_fix_dict['V']]
        
        # 将固定空位索引添加到数组
        vacancy_ind = np.concatenate([vacancy_ind, fix_ind_dict['V']]).astype(int)

        # 从骨架中删除空位
        vacancy_bool = np.array([True] * len(backbone))
        vacancy_bool[vacancy_ind] = False
        backbone = backbone[vacancy_bool]

        # 更新固定原子字典
        del fix_ind_dict['V']
        for key, value in fix_ind_dict.items():
            for i in range(len(value)):
                value[i] -= np.sum(vacancy_ind < value[i])
            fix_ind_dict[key] = value
    
    fix_ind_arr = np.concatenate(list(fix_ind_dict.values())).astype(np.int32)

    # 初始化新原子位置
    new_posi = []
    len_atoms = len(backbone)
    old_posi = backbone.get_positions()
    cell = backbone.get_cell()
    cell_inv = np.linalg.inv(cell)

    while len(new_posi) < len_atoms:
        # 记录循环开始时间
        time_i = time()
        
        while True:
            # 附加一个新的原子
            new_posi.append(old_posi[len(new_posi)].copy())
            
            # 给新原子位置添加偏移
            if not pin_the_fixed or len(new_posi) - 1 not in fix_ind_arr:
                direc_vec = np.random.rand(3) - 0.5
                direc_vec /= np.linalg.norm(direc_vec)
                new_posi[-1] += direc_vec * ran_radi

            # 获取最新原子的最小距离
            rel_new_posi = np.array(new_posi) @ cell_inv
            if len(new_posi) != 1:
                min_dist = np.min(np.linalg.norm(((rel_new_posi[:-1] - rel_new_posi[-1] + np.array([0.5]*3)) % 1.0 - np.array([0.5]*3)) @ cell, axis=1))
            else:
                min_dist = cutoff_r + 1.

            # 获取循环已用时间
            time_f = time()
            time_d = time_f - time_i

            # 如果最新原子正确定位，则跳出循环
            if min_dist > cutoff_r:
                if log and len(new_posi) % 100 == 0:
                    print(f"( {len(new_posi)} th / {len_atoms} ) new atom position found")
                break

            # 如果最新原子距离过近，则移除最新原子
            elif time_d < max_trial_sec:
                new_posi.pop()

            # 如果尝试超时，则从头开始
            else:
                new_posi = []
                break
    
    # 创建新原子对象
    new_atoms = backbone.copy()
    new_atoms.set_positions(new_posi, apply_constraint=False)

    # 打乱原子位置
    shuffle_ind = np.setdiff1d(range(len(new_atoms)), np.concatenate(list(fix_ind_dict.values())), True)
    new_positions = new_atoms.get_positions().copy()
    new_positions[shuffle_ind] = np.random.permutation(new_positions[shuffle_ind])
    new_atoms.set_positions(new_positions, apply_constraint=False)
    
    # 校正化学符号
    new_species = np.array(['XX'] * len(new_atoms))
    new_species[shuffle_ind] = shffl_spec_list
    if len(fix_ind_arr):
        new_species[fix_ind_arr] = count2list(num_fix_dict)

    # 设置化学符号
    new_atoms.set_chemical_symbols(new_species)
    
    # 按化学序号排序
    if sort_atoms:
        new_atoms = new_atoms[np.argsort(new_atoms.get_atomic_numbers())]

    return new_atoms
