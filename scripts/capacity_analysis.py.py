#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
用法: python capacity_analysis.py <路径> <阳离子1> <阳离子2> ... <阳离子n> <阳离子数量>

该脚本从指定路径中读取结构文件，计算给定阳离子的最大重量比容量（mAh/g），并将结果保存为capacity.csv。
"""

import os
import pandas as pd
from analyzer import BatteryAnalyzer
from pymatgen.io.vasp.inputs import Structure
from sys import argv

def calculate_capacity(path, cations, num_cations):
    """
    计算最大重量比容量，并将结果保存为CSV文件。

    参数:
        path: 结构文件的路径
        cations: 阳离子列表
        num_cations: 阳离子数量
    """
    slab_files = [os.path.join(path, i) for i in os.listdir(path)]
    cap_dict = {}

    for ion in cations:
        ion_dict = {}
        for file in slab_files:
            structure = Structure.from_file(file)
            slab_name = structure.formula.replace(' ', '')
            ele_list = structure.composition.elements

            oxi_dict = {ele.name: +1 for ele in ele_list}
            oxi_dict["B"] = -2

            structure.add_oxidation_state_by_element(oxi_dict)
            battery = BatteryAnalyzer(struc_oxid=structure, cation=ion)
            capacity = battery.get_max_capgrav(num_cations)  # 返回最大重量比容量 (mAh/g)
            ion_dict[slab_name] = round(capacity, 2)
        cap_dict[ion] = ion_dict

    df = pd.DataFrame.from_dict(cap_dict, orient='index').T
    df.to_csv('capacity.csv')

if __name__ == "__main__":
    if len(argv) < 4:
        print("用法: python capacity_analysis.py <路径> <阳离子1> <阳离子2> ... <阳离子n> <阳离子数量>")
    else:
        path = argv[1]
        cations = argv[2:-1]
        num_cations = int(argv[-1])
        calculate_capacity(path, cations, num_cations)
