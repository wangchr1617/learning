"""
用法: python vasp2xyz.py [label]

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
    if os.path.exists('screen_tmp'):
        os.remove('screen_tmp')

    # 删除NEP-dataset.xyz文件（如果存在）
    if os.path.exists('NEP-dataset.xyz'):
        os.remove('NEP-dataset.xyz')
    
    # 逐行读取xmllist中的每个文件路径
    with open('xmllist', 'r') as file_list:
        for line in file_list:
            xml = line.strip('\n')
            print(f"Processing file: {xml}")

            try:
                # 读取vasprun.xml中的原子结构
                b = read(xml, index=":")
            except Exception:
                # 如果读取vasprun.xml失败，则尝试读取OUTCAR文件
                outcar_path = xml.replace("vasprun.xml", "OUTCAR")
                b = read(outcar_path, index=":")
                print(f"Fallback to: {outcar_path}")
        
            # 检查每个离子步骤的收敛性
            os.system(f"grep -B 1 E0 {xml.replace('vasprun.xml','OSZICAR')} | grep -E 'DAV|RMM' | awk '{{if($2>=120) print 0; else print 1}}' > screen_tmp")
            screen = np.loadtxt("screen_tmp")
        
            # 如果screen是标量，将其转换为列表
            if screen.ndim == 0:
                screen = [screen]
        
            # 处理每个收敛的离子步骤
            for ind, is_converged in enumerate(screen):
                if is_converged == 1:
                    # 计算应力张量和体积
                    xx, yy, zz, yz, xz, xy = -b[ind].calc.results['stress'] * b[ind].get_volume()
                    b[ind].info['virial'] = np.array([(xx, xy, xz), (xy, yy, yz), (xz, yz, zz)])
                    
                    # 删除应力信息
                    del b[ind].calc.results['stress']
                    b[ind].pbc = True
                    
                    # 添加配置类型标签
                    b[ind].info['config_type'] = label

                    # 将原子结构写入XYZ文件
                    write("NEP-dataset.xyz", b[ind], append=True)
    
    # 删除临时文件
    os.remove('screen_tmp')
    os.remove('xmllist')

if __name__ == "__main__":
    # 检查命令行参数数量
    if len(sys.argv) == 2:
        label = sys.argv[1]  # 使用提供的标签
    else:
        label = 'low'  # 默认标签
    main(label)
