# Usage: python phonbycalorine.py nep|dft|both
# Necessary files in ./: POSCAR PRIMCELL.vasp nep.txt
import numpy as np
import os
import pandas as pd
import phonopy
import sys
from ase import Atoms
from ase.calculators.mixing import SumCalculator
from ase.calculators.singlepoint import SinglePointCalculator
from ase.io import read
from calorine.calculators import CPUNEP
from calorine.tools import relax_structure
from dftd3.ase import DFTD3
from matplotlib import pyplot as plt
from phonopy import Phonopy
from phonopy.structure.atoms import PhonopyAtoms
from phonopy.units import THzToCm
from seekpath import get_explicit_k_path
from typing import Any, Dict
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'

if len(sys.argv) < 2:
    print("Usage: python phonbycalorine.py [nep|dft|both]")
    sys.exit(1)

def cal_prim_matrix(unit_cell, prim_cell):
    unit_inv = np.linalg.inv(unit_cell)
    pc_matrix = np.dot(unit_inv, prim_cell)
    return pc_matrix

def NEP_calculator(dftd3=False):
    '''
    需要考虑色散校正时设置dftd3=True
    '''
    if os.path.exists("./nep.txt"):
        if dftd3 == True:
            calculator = SumCalculator([CPUNEP("./nep.txt"), DFTD3(method="pbe", damping="d3bj")]) # 默认是DFTD3-BJ
        else:
            calculator = CPUNEP("./nep.txt")
        return calculator
    else:
        print("The file nep.txt does not exist.")
        sys.exit(1)
        
def forcesbynep(structure):    
    structure.calc = NEP_calculator()
    forces = structure.get_forces().copy()
    return forces
    
def fcbynep(structure, pc_matrix, sc_matrix):    
    structure.calc = NEP_calculator()
    # relax_structure(structure, fmax=0.0001) # 默认不再弛豫结构
    structure_ph = PhonopyAtoms(symbols=structure.symbols, cell=structure.cell, scaled_positions=structure.get_scaled_positions(), pbc=structure.pbc)
    phonon = Phonopy(structure_ph, sc_matrix, pc_matrix)
    phonon.generate_displacements(distance=0.01) # 默认位移距离是0.01埃.
    forces = []
    for structure_ph in phonon.supercells_with_displacements:
        structure_ase = Atoms(symbols=structure_ph.symbols, cell=structure_ph.cell, scaled_positions=structure_ph.scaled_positions, pbc=structure.pbc)
        f = forcesbynep(structure_ase)
        forces.append(f)
    phonon.forces = forces
    phonon.produce_force_constants()
    # phonon.set_force_constants_zero_with_radius(cutoff_radius=20) # 设两倍NEP势函数截断半径外的力常数为0
    # phonon.save(settings={'force_constants': True}) # 用于保存phonopy信息
    return phonon
    
def fcbydft(structure, pc_matrix, sc_matrix):
    if os.path.exists("phonopy.yaml"):
        phonon = phonopy.load("phonopy.yaml") # 通过phonopy.yaml加载DFT+phonopy的计算结果
    else:
        if os.path.exists("FORCE_SETS") or os.path.exists("FORCE_CONSTANTS"):
            phonon = phonopy.load(supercell_matrix=sc_matrix, 
                                  primitive_matrix=pc_matrix,
                                  unitcell_filename="POSCAR") # 通过FORCE_SETS或FORCE_CONSTANTS加载DFT+phonopy的计算结果
        else:
            print("The file FORCE_SETS or FORCE_CONSTANTS does not exist.")
            return None
    return phonon
    
def phonDataFrame(structure, phonon):
    structure_tuple = (structure.cell, structure.get_scaled_positions(), structure.numbers)
    path = get_explicit_k_path(structure_tuple) # 自动查找高对称点
    labels = path['explicit_kpoints_labels']
    labels = ['$\Gamma$' if m == 'GAMMA' else m for m in labels]
    labels = [m.replace('_', '$_') + '$' if '_' in m else m for m in labels]
    phonon.run_band_structure([path['explicit_kpoints_rel']], with_eigenvectors=True, with_group_velocities=True) # 额外输出本征值和群速度
    band = phonon.get_band_structure_dict() # dict_keys = ['qpoints', 'distances', 'frequencies', 'eigenvectors', 'group_velocities']
    df = pd.DataFrame(band['frequencies'][0])
    df.index = path['explicit_kpoints_linearcoord']
    return df, path, labels
    
def plt_phondisp(df1, g1, path, labels, figname, df2=None, g2=None):
    fig, ax = plt.subplots(figsize=(4.2, 3), dpi=140)
    for i, col in enumerate(df1.columns):
        if i == 0:
            ax.plot(df1.index, df1[col], c="orange", label=g1)
        else:
            ax.plot(df1.index, df1[col], c="orange")
    if df2 is not None:
        for i, col in enumerate(df2.columns):
            if i == 0:
                ax.plot(df2.index, df2[col], c="purple", alpha=0.5, ls="--", lw=0.8, label=g2)
            else:
                ax.plot(df2.index, df2[col], c="purple", alpha=0.5, ls="--", lw=0.8)
    ax.set_xlim(df1.index.min(), df1.index.max())
    ax.set_xlabel('KPATH')
    ax.set_ylabel('Frequency (THz)')
    ax2 = ax.twinx()
    ax2.set_ylabel('Frequency (cm$^{-1}$)')
    ax2.set_ylim(THzToCm * np.array(ax.get_ylim()))
    df_path = pd.DataFrame(dict(labels=labels, positions=path['explicit_kpoints_linearcoord']))
    df_path.drop(df_path.index[df_path.labels == ''], axis=0, inplace=True)
    ax.set_xticks(df_path.positions)
    ax.set_xticklabels(df_path.labels)
    for xp in df_path.positions:
        ax.axvline(xp, c='black', alpha=0.5, ls="--", lw=0.8)
    ax.axhline(y=0.0, c='black', alpha=0.5, ls="--", lw=0.8)
    ax.legend()
    plt.tight_layout()
    plt.savefig(figname, bbox_inches='tight')

def cal_dos(phonon):
    phonon.run_mesh([30, 30, 30], with_eigenvectors=False, with_group_velocities=False, is_gamma_center=True) 
    phonon.run_total_dos(freq_pitch=0.02)
    phonon.run_thermal_properties(t_min=0, t_max=500, t_step=5)
    dos = phonon.get_total_dos_dict()
    tprop = phonon.get_thermal_properties_dict()
    return dos, tprop
    
def plt_phondos(p1, g1, figname, p2=None, g2=None):    
    fig, ax = plt.subplots(figsize=(4.2, 3), dpi=140)
    dos1, _ = cal_dos(p1)
    ax.plot(dos1['frequency_points'], dos1['total_dos'], c="orange", label=g1)
    if p2 is not None:
        dos2 = cal_dos(p2)
        ax.plot(dos2['frequency_points'], dos2['total_dos'], c="purple", alpha=0.5, ls="--", label=g2)
    ax.set_xlabel('Frequency (THz)')
    ax.set_ylabel('Density of states')
    plt.tight_layout()
    plt.legend()
    plt.savefig(figname, bbox_inches='tight')

def plt_phontprop():
    fig, ax = plt.subplots(figsize=(4.2, 3), dpi=140)
    _, tprop = cal_dos(p1)
    ax.plot(tprop['temperatures'], tprop['heat_capacity'], c="orange", label="Heat_capacity")
    ax.set_xlabel('Temperature')
    ax.set_ylabel('Heat capacity')
    ax.set_xlim(0,)
    plt.tight_layout()
    plt.legend()
    plt.savefig("tprop.png", bbox_inches='tight')

fctype = sys.argv[1] # 判断任务类型
if fctype not in ["nep", "dft", "both"]:
    print("Error: The fctype must be 'nep', 'dft' or 'both'.")
    sys.exit(1)
unit_cell = read("./POSCAR")
# prim_cell = read("./PRIMCELL.vasp") # 使用VASPKIT 602查找原胞
# pc_matrix = cal_prim_matrix(unit_cell.get_cell(), prim_cell.get_cell()) # 根据超胞和原胞计算原胞矩阵
pc_matrix = 'auto' # 一般用auto即可
sc_matrix = [4, 4, 4] # 超胞矩阵
try:
    g2 = None
    p2 = None
    if fctype == "nep":
        g1 = "nep"
        p1 = fcbynep(unit_cell, pc_matrix, sc_matrix)
    elif fctype == "dft":
        g1 = "dft"
        p1 = fcbydft(unit_cell, pc_matrix, sc_matrix)
    elif fctype == "both":
        g1 = "nep"
        g2 = "dft"
        p1 = fcbynep(unit_cell, pc_matrix, sc_matrix)
        p2 = fcbydft(unit_cell, pc_matrix, sc_matrix)
        
    if p1 is None:
        print(f"Error: Phonon object not created. Check the fctype ('{fctype}') or input files.") 
        sys.exit(1)
    else:
        df1, path, labels = phonDataFrame(unit_cell, p1)
        if p2 is not None:
            df2, path, labels = phonDataFrame(unit_cell, p2)
            plt_phondisp(df1, g1, path, labels, f"{fctype}_disp.png", df2=df2, g2=g2)
            plt_phondos(p1, g1, f"{fctype}_dos.png", p2=p2, g2=g2)
        else:
            plt_phondisp(df1, g1, path, labels, f"{fctype}_disp.png")
            plt_phondos(p1, g1, f"{fctype}_dos.png")     
            plt_phontprop()
except Exception as e:
    print(f"An error occurred: {e}")
    sys.exit(1)
