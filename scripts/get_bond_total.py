# Usage: python get_bond_total.py
import os
import sys
from pprint import pprint
from collections import defaultdict
from itertools import combinations, product

import numpy as np
from pymatgen.core.structure import Structure

info="""使用说明:
python get_bond_total.py
不需要输入任何参数 
"""
print(info)

struct = Structure.from_file("POSCAR")

# 获得指定的两组元素的坐标
sites_custom1 = struct.sites
index_custom1 = [idx+1 for idx, site in enumerate(struct.sites)] # 因为python的索引编号是从0开始的

pairs = list(combinations(sites_custom1, r=2))
idxs = list(combinations(index_custom1, r=2))
# 获得每个原子对的距离并且存储
pairs_dist = defaultdict(list)
for pair, idx in zip(pairs, idxs):
    d = struct.lattice.get_all_distances(pair[0].frac_coords, pair[1].frac_coords)[0][0]
    pair_name = str(pair[0].specie) + str(idx[0]) + '-' + str(pair[1].specie) + str(idx[1])
    pairs_dist[np.round(d, decimals=3)].append(pair_name)
pairs_dist = sorted(pairs_dist.items())

# 打印相关信息
print("Note: --------------------")
print("{} {}".format("Number", "Elements"))
for idx, site in enumerate(struct.sites):
    print("{:<8} {:<8}".format(idx+1, site.specie))

print("Note: --------------------")
print("{} {} {}".format("Species", "Species", "distance"))
for d, pairs in pairs_dist:
    print("distance = {}, Number={}".format(d, len(pairs)))
    print("    "+"   ".join(pairs))

with open("distance.dat", "w") as f:
    for d, pairs in pairs_dist:
        f.write("distance = {}, Number={}\n".format(d, len(pairs)))
        f.write("    "+"   ".join(pairs))
        f.write("\n")
