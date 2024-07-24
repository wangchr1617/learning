#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
用法: python vasp_to_xyz.py [label]

该脚本读取VASP输出文件（vasprun.xml或OUTCAR）并将其转换为XYZ格式。
通过使用OSZICAR筛选收敛的离子步骤，并且仅处理收敛的离子步骤。
结果保存在NEP-dataset.xyz文件中，并附加了关于应力张量和配置类型的附加信息。
"""

from ase.io import read, write
import numpy as np
import sys
import os

def main(label='low'):
    # 查找当前目录及其子目录下的vasprun.xml文件，并将文件路径保存到xmllist
    os.system("find . -name vasprun.xml > xmllist")
    # 删除临时文件screen_tmp（如果存在）
    os.system("if [ -f 'screen_tmp' ]; then rm screen_tmp; fi")
    # 删除NEP-dataset.xyz文件（如果存在）
    os.system("if [ -f 'NEP-dataset.xyz' ]; then rm NEP-dataset.xyz; fi")
    
    # 逐行读取xmllist中的每个文件路径
    for line in open('xmllist'):
        xml = line.strip('\n')
        print(xml)
        try:
            b = read(xml, index=":")
        except Exception:
            # 如果读取vasprun.xml失败，则尝试读取OUTCAR文件
            b = read(xml.replace("vasprun.xml", "OUTCAR"), index=":")
            print(xml.replace("vasprun.xml", "OUTCAR"))
        
        # 检查每个离子步骤的收敛性
        os.system(f"grep -B 1 E0 {xml.replace('vasprun.xml','OSZICAR')} | grep -E 'DAV|RMM' | awk '{{if($2>=120) print 0; else print 1}}' > screen_tmp")
        screen = np.loadtxt("screen_tmp")
        
        # 如果screen是标量，将其转换为列表
        if screen.ndim == 0:
            screen = [screen]
        
        # 处理每个收敛的离子步骤
        for ind, i in enumerate(screen):
            if i == 1:
                xx, yy, zz, yz, xz, xy = -b[ind].calc.results['stress'] * b[ind].get_volume()
                b[ind].info['virial'] = np.array([(xx, xy, xz), (xy, yy, yz), (xz, yz, zz)])
                del b[ind].calc.results['stress']
                b[ind].pbc = True
                b[ind].info['config_type'] = label
                write("NEP-dataset.xyz", b[ind], append=True)
    
    # 删除临时文件
    os.system("rm screen_tmp")
    os.system("rm xmllist")

if __name__ == "__main__":
    if len(sys.argv) == 2:
        label = sys.argv[1]
    else:
        label = 'low'
    main(label)
