"""
用法: python calculate_rdf.py

功能描述:
该脚本用于计算指定区域内 A-B 原子对的径向分布函数 (RDF)。
用户可以通过选择不同的区域来动态调整所计算的原子对。

参数说明:
- 输入文件:
  - 'npt0.tpr': 拓扑文件。
  - '10-20ns.xtc': 轨迹文件。

结果:
- `rdf.txt`: 包含 RDF 结果的文本文件。

依赖库:
- numpy
- MDAnalysis
"""

import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis.rdf import InterRDF

def main():
    # 创建一个MDAnalysis的Universe对象，包含拓扑和轨迹
    u = mda.Universe('npt0.tpr', '10-20ns.xtc')

    # 选择A和B原子，选择条件为z坐标在30到36之间
    select_A = u.select_atoms('name OW and (prop z > 30 and prop z < 36)', updating=True)
    select_B = u.select_atoms('name OW and (prop z > 30 and prop z < 36)', updating=True)

    # 计算A-B原子的RDF
    rdf = InterRDF(select_A, select_B, range=(0.0, 20.0), nbins=200, n_threads=12)
    rdf.run()

    # 设置第一个RDF值为0，避免数值不稳定
    rdf.results.rdf[0] = 0.0

    # 将RDF结果保存为文本文件
    np.savetxt("rdf.txt", np.column_stack((rdf.results.bins, rdf.results.rdf)), header="r (Angstrom)\t g(r)")
    print("RDF计算完成，结果已保存至 rdf.txt 文件中。")

if __name__ == "__main__":
    main()
