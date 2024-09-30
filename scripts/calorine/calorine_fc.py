import numpy as np
import os
import shutil
import sys
from ase import Atoms
from ase.calculators.mixing import SumCalculator
from ase.io import read, write
from calorine.calculators import CPUNEP
from calorine.tools import relax_structure
from dftd3.ase import DFTD3
from phono3py import Phono3py
from phono3py.file_IO import write_fc2_to_hdf5, write_fc3_to_hdf5
from phono3py.interface.phono3py_yaml import Phono3pyYaml
from phonopy.structure.atoms import PhonopyAtoms

def NEP_calculator(dftd3=False):
    if os.path.exists("./nep.txt"):
        if dftd3:
            calculator = SumCalculator([CPUNEP("./nep.txt"), DFTD3(method="pbe", damping="d3bj")])  # 默认是DFTD3-BJ
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

def fcbynep(structure, pc_matrix, sc_matrix, ps_matrix):
    structure.calc = NEP_calculator()
    relax_structure(structure, fmax=1e-5)  # 弛豫结构
    write('POSCAR', structure)
    structure_ph = PhonopyAtoms(symbols=structure.symbols, cell=structure.cell, scaled_positions=structure.get_scaled_positions(), pbc=structure.pbc)
    phonon = Phono3py(structure_ph, supercell_matrix=sc_matrix, primitive_matrix=pc_matrix, phonon_supercell_matrix=ps_matrix)
    phonon.generate_displacements(distance=0.03, cutoff_pair_distance=None, is_plusminus="auto", is_diagonal=True)
    # phonon.generate_fc2_displacements(distance=0.01, is_plusminus="auto", is_diagonal=False)
    phonon.save("phono3py_disp.yaml")
    
    # 计算二阶力常数对应的超胞结构受力
    print("Number of displacements for special fc2: ", len(phonon.phonon_supercells_with_displacements))
    fc2 = []
    it = 1
    for structure_ph in phonon.phonon_supercells_with_displacements:
        structure_ase = Atoms(symbols=structure_ph.symbols, cell=structure_ph.cell, scaled_positions=structure_ph.scaled_positions, pbc=structure.pbc) 
        # write(f"POSCAR_FC2-{it:05d}", structure_ase)
        f = forcesbynep(structure_ase)
        fc2.append(f)
        it += 1
    phonon.phonon_forces = fc2
    fc2 = np.array(fc2).reshape(-1, 3)
    np.savetxt('FORCES_FC2', fc2)
    
    # 计算三阶力常数对应的超胞结构受力
    print("Number of displacements: ", len(phonon.supercells_with_displacements))
    fc3 = []
    it = 1
    for structure_ph in phonon.supercells_with_displacements:
        structure_ase = Atoms(symbols=structure_ph.symbols, cell=structure_ph.cell, scaled_positions=structure_ph.scaled_positions, pbc=structure.pbc)
        # write(f"POSCAR-{it:05d}", structure_ase)
        f = forcesbynep(structure_ase)
        fc3.append(f)
        it += 1
    phonon.forces = fc3
    fc3 = np.array(fc3).reshape(-1, 3)
    np.savetxt('FORCES_FC3', fc3)
    
    force_constants_2 = phonon.produce_fc2()
    write_fc2_to_hdf5(force_constants_2, filename="fc2.hdf5", p2s_map=None, physical_unit=None, compression="gzip")
    force_constants_3 = phonon.produce_fc3()
    write_fc3_to_hdf5(force_constants_3, filename="fc3.hdf5", p2s_map=None, compression="gzip")
    
root_dir = os.getcwd()  # here we set the root directory
work_dir = os.path.join(root_dir, 'phono3py_workdir')
shutil.rmtree(work_dir, ignore_errors=True)  # Deletes current working dir
os.mkdir(work_dir)
os.chdir(work_dir)
shutil.copyfile('../nep.txt', 'nep.txt')

unit_cell = read(os.path.join(root_dir, 'POSCAR'))
pc_matrix = 'auto'
sc_matrix = [4, 4, 2]  # fc3可以使用比fc2更小的超胞
ps_matrix = [4, 4, 2]  # fc2
p1 = fcbynep(unit_cell, pc_matrix, sc_matrix, ps_matrix)
print("Over!!!")

