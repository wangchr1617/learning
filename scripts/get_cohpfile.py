import os
import sys
from pprint import pprint
from collections import defaultdict
from itertools import combinations, product

import numpy as np
from pymatgen.core.structure import Structure

info="""使用说明:
python get_pairs.py [COHPstartEnergy] [COHPendEnergy] [element1] [element2] [lower_d_limit] [upper_d_limit]
指定其实计算的能量COHPstartEnergy和终止计算的能量COHPendEnergy
指定两种元素element1和element2, 两种元素可以相同, 效果是求两种元素在单胞内所有原子间距离
指定原子间距离，获得指定距离下原子对。
"""
print(info)

COHPstartEnergy = sys.argv[1]
COHPendEnergy   = sys.argv[2]
species_custom1 = sys.argv[3]
species_custom2 = sys.argv[4]

try:
    lower_d = float(sys.argv[5])
    upper_d = float(sys.argv[6])
except:
    lower_d = 0.0
    upper_d = 5.0


struct = Structure.from_file("POSCAR")

with open("lobsterinIsmear_5", "w") as f:
    f.write('COHPstartEnergy  {}\n'.format(COHPstartEnergy)) 
    f.write('COHPendEnergy    {}\n'.format(COHPendEnergy))
    f.write('usebasisset pbeVaspFit2015\n') # 基组（直接使用根据vasp拟合的基组以及默认的基函数）
    for spe in struct.types_of_specie:
        f.write('basisfunctions {}\n'.format(spe.name)) # 基组（直接使用根据vasp拟合的基组以及默认的基函数）
    # 这里会出现一个非常严重的计算问题！！！！！！！！！！！！！！
    # lobster认为你指定的原子对是不具有周期性的
    # 你用pymatgen脚本找到的距离是包含周期性的，把这原子对输入给lobsterin
    # 它认不出来这个距离是周期性的，它会按照原胞内的距离考虑两个原子的成键。
    # 所以这里我抛弃了设置原子对来计算成键强度的方法。
    # 改用设置键长来获得原子对，lobster有自己的算法来获得原子对。
    # for pair, idx, d in pairs_idxs_d:
    #     if lower_d <= d <= upper_d:
    #         f.write("cohpbetween atom {} and atom {}\n".format(idx[0], idx[1]))
    f.write("cohpGenerator from {} to {} type {} type {}\n".format(lower_d, upper_d, species_custom1, species_custom2))

with open("lobsterinIsmear_0", "w") as f:
    f.write('COHPstartEnergy  {}\n'.format(COHPstartEnergy)) 
    f.write('COHPendEnergy    {}\n'.format(COHPendEnergy))
    f.write('usebasisset pbeVaspFit2015\n') # 基组（直接使用根据vasp拟合的基组以及默认的基函数）
    f.write('gaussianSmearingWidth 0.05\n')
    for spe in struct.types_of_specie:
        f.write('basisfunctions {}\n'.format(spe.name)) # 基组（直接使用根据vasp拟合的基组以及默认的基函数）
    # 这里会出现一个非常严重的计算问题！！！！！！！！！！！！！！
    # lobster认为你指定的原子对是不具有周期性的
    # 你用pymatgen脚本找到的距离是包含周期性的，把这原子对输入给lobsterin
    # 它认不出来这个距离是周期性的，它会按照原胞内的距离考虑两个原子的成键。
    # 所以这里我抛弃了设置原子对来计算成键强度的方法。
    # 改用设置键长来获得原子对，lobster有自己的算法来获得原子对。
    # for pair, idx, d in pairs_idxs_d:
    #     if lower_d <= d <= upper_d:
    #         f.write("cohpbetween atom {} and atom {}\n".format(idx[0], idx[1]))
    f.write("cohpGenerator from {} to {} type {} type {}\n".format(lower_d, upper_d, species_custom1, species_custom2))



dirs = "{}_{}_{}_{}".format(species_custom1, species_custom2, lower_d, upper_d)
if not os.path.exists(dirs):
    os.mkdir(dirs)
print("Note: ------------------------------------------------------")
print("    You need to do it before running Lobster-4.1.0")
print("    cp lobsterraw lobsterin")
print("    You need to do it after running Lobster-4.1.0")
print("    cp {COBICAR.lobster,COHPCAR.lobster,COOPCAR.lobster,ICOBILIST.lobster,ICOHPLIST.lobster,ICOOPLIST.lobster,lobsterin}" + "  {}".format(dirs))
print("------------------------------------------------------------")
