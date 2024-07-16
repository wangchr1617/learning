import numpy as np
import os
from ase.build import sort
from ase.io import read, write
from RAG import random_atoms_gen as rag

def rag_iter(structure,
             idx,
             num_spec_dict, 
             fix_ind_dict,
             iter_num=20, 
             random_radi=0.01,
             cutoff_frac=0.9,
             strain_max=1.03, 
             strain_min=0.97, 
             vacuum_max=0.0):
    for i in range(iter_num):
        strain_ratio = np.random.rand(3) * (strain_max - strain_min) + strain_min
        vacuum = [0., 0., np.random.rand() * vacuum_max]
        atoms = rag(structure,
                    num_spec_dict = num_spec_dict,
                    fix_ind_dict  = fix_ind_dict,
                    cutoff_frac   = cutoff_frac,
                    random_radi   = random_radi,
                    strain_ratio  = strain_ratio,
                    vacuum        = vacuum,
                    log           = True
        )
        write('./rag_structures/{}_{}-{:0>2d}.vasp'.format(idx, random_radi, i), sort(atoms))

rag_path = './rag_structures'
if not os.path.exists(rag_path):
    os.makedirs(rag_path)

backbone_path = './backbone/'
radi_list = [0.08, 0.16, 0.32, 0.64]
idx = '6'
backbone = backbone_path + '6.vasp'
structure = read(backbone)
num_spec_dict = {'Ge':30, 'Te':32}
fix_ind = np.arange(30, 62)
for i in range(len(radi_list)):
    fix_ind_dict = {'Te': fix_ind}
    rag_iter(structure, idx, num_spec_dict, fix_ind_dict, iter_num=20, random_radi=radi_list[i])

