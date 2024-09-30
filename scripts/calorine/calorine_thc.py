import os
import shutil
from glob import glob

import numpy as np
import h5py
from ase.io import read, write
from calorine.calculators import CPUNEP
from calorine.tools import relax_structure
from dftd3.ase import DFTD3
from matplotlib import pyplot as plt

def NEP_calculator(dftd3=False):
    """初始化 NEP 计算器，可以选择性地添加 DFTD3。"""
    if os.path.exists("./nep.txt"):
        if dftd3:
            calculator = SumCalculator([CPUNEP("./nep.txt"), DFTD3(method="pbe", damping="d3bj")])
        else:
            calculator = CPUNEP("./nep.txt")
        return calculator
    else:
        print("The file nep.txt does not exist.")
        sys.exit(1)

def main(input_file, dim, mesh, T_min, T_max, T_step, work_dir_name='phono3py_workdir'):
    # 设置目录
    root_dir = os.getcwd()
    work_dir = os.path.join(root_dir, work_dir_name)
    shutil.rmtree(work_dir, ignore_errors=True)  # 删除已有的工作目录
    os.mkdir(work_dir)
    os.chdir(work_dir)
    shutil.copyfile(os.path.join(root_dir, 'nep.txt'), 'nep.txt')

    print(f'Root directory: {root_dir}')
    print(f'Working directory: {work_dir}')

    # 读取结构并优化
    prim = read(os.path.join(root_dir, input_file))
    prim.calc = NEP_calculator()
    relax_structure(prim, fmax=1e-5)

    # 生成位移用于力常数计算
    cmd = f'phono3py -d --dim="{dim[0]} {dim[1]} {dim[2]}" --dim-fc2="{dim[0]} {dim[1]} {dim[2]}"'
    write('POSCAR', prim)
    print(f'Running command: {cmd}')
    os.system(cmd)

    # 计算 FC2 的力
    fnames = sorted(glob('POSCAR_FC2-*'))
    forces_data = []
    for it, fname in enumerate(fnames):
        structure = read(fname)
        structure.calc = NEP_calculator()
        forces = structure.get_forces()
        forces_data.append(forces)
        print(f'FC2: Calculating supercell {it:5d} / {len(fnames)}, f_max {np.max(np.abs(forces)):8.5f}')
    forces_data = np.array(forces_data).reshape(-1, 3)
    np.savetxt('FORCES_FC2', forces_data)

    # 计算 FC3 的力
    fnames = sorted(glob('POSCAR-*'))
    forces_data = []
    for it, fname in enumerate(fnames):
        structure = read(fname)
        structure.calc = NEP_calculator()
        forces = structure.get_forces()
        forces_data.append(forces)
        if it % 100 == 0:
            print(f'FC3: Calculating supercell {it:5d} / {len(fnames)}, f_max= {np.max(np.abs(forces)):8.5f}')
    forces_data = np.array(forces_data).reshape(-1, 3)
    np.savetxt('FORCES_FC3', forces_data)

    # 计算力常数
    cmd = f'phono3py --cfc --hdf5-compression gzip --dim="{dim[0]} {dim[1]} {dim[2]}" --dim-fc2="{dim[0]} {dim[1]} {dim[2]}"'
    print(f'Running command: {cmd}')
    os.system(cmd)

    # 计算热导率
    cmd = (f'phono3py -q --fc2 --fc3 --bterta --dim="{dim[0]} {dim[1]} {dim[2]}" '
           f'--dim-fc2="{dim[0]} {dim[1]} {dim[2]}" --tmin={T_min} --tmax={T_max} '
           f'--tstep={T_step} --mesh "{mesh[0]} {mesh[1]} {mesh[2]}"')
    print(f'Running command: {cmd}')
    os.system(cmd)
    os.chdir(root_dir)

    # 读取结果
    with h5py.File(f'{work_dir}/kappa-m{mesh[0]}{mesh[1]}{mesh[2]}.hdf5') as fobj:
        temperatures = fobj['temperature'][:]
        kappa = fobj['kappa'][:]
        gamma = fobj['gamma'][:]
        freqs = fobj['frequency'][:]

    # 绘制热导率
    fig, ax = plt.subplots(figsize=(3.8, 3), dpi=140)
    ax.plot(temperatures, kappa[:, 0], '-', label=r'$\kappa_{xx}$')
    ax.plot(temperatures, kappa[:, 1], '--', label=r'$\kappa_{yy}$')
    ax.plot(temperatures, kappa[:, 2], '-', label=r'$\kappa_{zz}$')
    ax.legend(loc='upper right', frameon=False)
    ax.set_xlim([50, temperatures.max()])
    ax.set_ylim(0, max([kappa[:, 0].max(), kappa[:, 1].max(), kappa[:, 2].max()]) * 1.1)
    ax.set_xlabel('Temperature (K)')
    ax.set_ylabel('Thermal conductivity (W/m/K)')
    fig.tight_layout()
    plt.savefig('T-T.png', bbox_inches='tight')

    # 绘制声子寿命
    fig, ax = plt.subplots(figsize=(3.8, 3), dpi=140)
    T_indices = [9, 29, 99]
    for T_index in T_indices:
        print(f'Temperature {temperatures[T_index]}')
        g = gamma[T_index].flatten()
        g_nonzero = np.where(g == 0, np.nan, g)
        lifetimes = np.where(g_nonzero > 0.0, 1.0 / (2 * 2 * np.pi * g_nonzero), np.nan)
        ax.semilogy(freqs.flatten(), lifetimes, 'o', label=f'T={temperatures[T_index]:.0f} K', alpha=0.5, markersize=2)

    ax.legend(loc='upper right', frameon=False)
    ax.set_xlabel('Phonon frequency (THz)')
    ax.set_ylabel('Phonon lifetime (ps)')
    ax.set_xlim([0, freqs.max() * 1.02])
    fig.tight_layout()
    plt.savefig('PL-PF.png', bbox_inches='tight')


if __name__ == "__main__":
    input_file = 'POSCAR'
    dim = (4, 4, 4)
    mesh = [16, 16, 16]
    T_min, T_max = 10, 1000
    T_step = 10
    main(input_file, dim, mesh, T_min, T_max, T_step)

