# 文件名: nose_hoover_mass_calculator.py
# 运行方法: python nose_hoover_mass_calculator.py -i POSCAR -T 300 -f 40 -u fs
# 功能描述: 根据给定的温度、频率等参数计算Nosé-Hoover链热浴质量（SMASS）

import argparse
import ase
import sys
import numpy as np
from ase.io import read


def nose_mass(temperature, ndof, t0, L):
    '''
    计算Nosé-Hoover链热浴质量（SMASS）。

    参数:
        temperature (float): 温度，以开尔文为单位。
        ndof (int): 自由度的数量。
        t0 (float): 振荡时间，以fs为单位。
        L (float): 第一基矢的长度，以Å为单位。

    返回:
        float: 计算的Nosé-Hoover链热浴质量。
    '''
    # Q in Energy * Time**2
    qtmp = (t0 * 1E-15 / np.pi / 2)**2 * \
        2 * ndof * ase.units.kB * temperature \
        * ase.units._e

    # Q in AMU * Å**2
    Q = qtmp / ase.units._amu / (L * 1E-10)**2

    return Q


def cnt_dof(atoms):
    '''
    计算系统的自由度数量。

    参数:
        atoms (ase.Atoms): ASE的Atoms对象，表示分子或晶体结构。

    返回:
        int: 自由度的数量。
    '''
    if atoms.constraints:
        from ase.constraints import FixAtoms, FixScaled, FixedPlane, FixedLine
        sflags = np.zeros((len(atoms), 3), dtype=bool)
        for constr in atoms.constraints:
            if isinstance(constr, FixScaled):
                sflags[constr.a] = constr.mask
            elif isinstance(constr, FixAtoms):
                sflags[constr.index] = [True, True, True]
            elif isinstance(constr, FixedPlane):
                mask = np.all(np.abs(np.cross(constr.dir, atoms.cell)) < 1e-5,
                              axis=1)
                if sum(mask) != 1:
                    raise RuntimeError(
                        'VASP要求FixedPlane约束的方向与一个晶胞轴平行')
                sflags[constr.a] = mask
            elif isinstance(constr, FixedLine):
                mask = np.all(np.abs(np.cross(constr.dir, atoms.cell)) < 1e-5,
                              axis=1)
                if sum(mask) != 1:
                    raise RuntimeError(
                        'VASP要求FixedLine约束的方向与一个晶胞轴平行')
                sflags[constr.a] = ~mask

        return np.sum(~sflags)
    else:
        return len(atoms) * 3 - 3


def parse_cml_args(cml):
    '''
    命令行参数解析器。

    参数:
        cml (list): 命令行参数列表。

    返回:
        argparse.Namespace: 解析后的命令行参数。
    '''
    arg = argparse.ArgumentParser(add_help=True)

    arg.add_argument('-i', dest='poscar', action='store', type=str,
                     default='POSCAR',
                     help='用于计算自由度数量的真实POSCAR文件')
    arg.add_argument('-u', dest='unit', action='store', type=str,
                     default='fs',
                     choices=['cm-1', 'fs'],
                     help='输入频率的默认单位')
    arg.add_argument('-T', '--temperature', dest='temperature',
                     action='store', type=float,
                     default=300,
                     help='温度，以开尔文为单位。')
    arg.add_argument('-f', '--frequency', dest='frequency',
                     action='store', type=float,
                     default=40,
                     help='温度振荡的频率')

    return arg.parse_args(cml)


if __name__ == '__main__':
    # 解析命令行参数
    arg = parse_cml_args(sys.argv[1:])

    # 将频率单位从cm-1转换为fs
    if arg.unit == 'cm-1':
        THzToCm = 33.3564095198152
        t0 = 1000 * THzToCm / arg.frequency
    else:
        t0 = arg.frequency

    # 读取POSCAR文件
    pos = read(arg.poscar)

    # 计算第一基矢的长度
    L = np.linalg.norm(pos.cell, axis=1)[0]

    # 计算自由度数量
    ndof = cnt_dof(pos)

    # 计算Nosé-Hoover链热浴质量
    Q = nose_mass(arg.temperature, ndof, t0, L)

    print("SMASS = {}".format(Q))
    