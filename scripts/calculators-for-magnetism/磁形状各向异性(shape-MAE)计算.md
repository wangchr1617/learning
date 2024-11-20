### 磁形状各向异性(shape-MAE)计算

磁各向异性(MA)是二维材料在有限温度下建立磁序的基础。磁各向异性能(total-MAE)分为自旋轨道耦合诱导的磁晶各向异性(SOC-MAE)，和磁偶极子诱导的形状各向异性(shape-MAE)。

磁偶极-偶极相互作用可通过下面公式来描述：

$$
E^{Dipole} = -\frac{1}{2}\frac{\mu_0}{4{\pi}}\sum_{i=1}^{N}\sum_{j=1}^{r_{max}}\frac{1}{r^3_{ij}}{\left[ \overrightarrow{M_1}\cdot\overrightarrow{M_2} - \frac{3}{r^2_{ij}}\left( \overrightarrow{M_1}\cdot\overrightarrow{r_{ij}} \right)\left( \overrightarrow{M_2}\cdot\overrightarrow{r_{ij}} \right) \right]}
$$

$\mu_0$ 是真空磁导率，$\overrightarrow{M_i}$ 是磁性原子在位点$i$的磁矩，$\overrightarrow{r_{ij}}$是位点$ i, j$ 的相对坐标矢量。

下面以单层CrSBr为例，计算形状各向异性

CrSBr的POSCAR如下：

```
CrSBr         
   1.00000000000000   
     3.5869562778239521    0.0000000000000000    0.0000000000000000
     0.0000000006263387    4.8226464297073965    0.0000000000000000
     0.0000000000000000    0.0000000000000000   25.6198997497999983
   Cr   S    Br
     2     2     2
Direct
  0.4999999939003753  0.9999999797715552  0.4594350922177504
  0.0000000059383842  0.5000000202991159  0.5405649192327502
  0.5000000060860624  0.4999999930217441  0.4767864028478482
  0.9999999939811991  0.0000000068330420  0.5232136124675983
  0.9999999915109647  0.9999999797482435  0.3889819115056300
  0.5000000085830142  0.5000000203263063  0.6110180987284154
```

CrSBr单胞有两个Cr原子，可以用ase读取POSCAR并通过ase.neighborlist获得近邻原子

```
import numpy as np
from datetime import datetime
from ase import Atoms
from ase.io import read
from ase.neighborlist import NeighborList

def read_cell(file_path, element, pbc):
    atom_object = read(file_path)
    atom_object.pbc = pbc
    element_list = atom_object.get_chemical_symbols()
    indices = [index for index, value in enumerate(element_list) if value != element]
    del atom_object[indices]
    return atom_object

def get_neighborlist(atom_object, cutoff):
    atom_num = len(atom_object)
    cutoff = [cutoff]*atom_num
    nl = NeighborList(cutoff, self_interaction=False, bothways=True)
    nl.update(atom_object)
    atom_neigh = []
    for i in range(atom_num):
        indices, offsets = nl.get_neighbors(i)
        dxyz = atom_object.positions[indices] + np.dot(offsets, atom_object.cell) - atom_object.positions[i]
        dr = np.sqrt(np.sum(dxyz**2, axis=1))
        atomR = np.column_stack((indices,dxyz,dr))
        atomR[:, 0] = atomR[:, 0].astype(int)
        atom_neigh.append(atomR)
    return atom_neigh

def get_mm(atom_object, magnetic_moment):
    atom_num = len(atom_object)
    mmx = [[0, 0, 0]]*atom_num
    mmy = [[0, 0, 0]]*atom_num
    mmz = [[0, 0, 0]]*atom_num
    for i in range(atom_num):
        mmx[i][0] = magnetic_moment[i]
        mmy[i][1] = magnetic_moment[i]
        mmz[i][2] = magnetic_moment[i]
    return mmx, mmy, mmz

def calculate_MSA(atom_neigh, magnetic_moment=[[0,0,1],[0,0,2]]):
    dipole_lst = []
  
    for i in range(len(atom_neigh)):
        mmi = magnetic_moment[i]
        dipole = 0
        for j in range(len(atom_neigh[i])):
            indices = atom_neigh[i][:,0][j].astype(int)
            mmj = magnetic_moment[indices]
            dxyz = atom_neigh[i][:,1:4][j]
            dr = atom_neigh[i][:,4][j]

            temp = np.dot(mmi, mmj)-3*np.dot(mmi, dxyz)*np.dot(mmj, dxyz)/(dr**2)
            temp_A = temp*((9.274*1e-24) ** 2)
            constant = 1e-7
            J_dipole = constant*temp_A/(2*dr**3)*1e30
            miueV_dipole = J_dipole * 6.241506363 * 1e24

            dipole = dipole + miueV_dipole
  
        dipole_lst.append(dipole)
    return dipole_lst
```

读取POSCAR，将截断距离 $r^{max}$设为1000Å，遍历范围内所有的近邻原子，分别计算磁矩沿x，y，z方向的磁偶极-偶极相互作用

```
if __name__ == '__main__':
    start_time = datetime.now()
    file_path = './POSCAR'
    pbc = [1, 1, 0]
    mag_element = 'Cr'
    mag_moment = [3, 3]
    rcut = 1000
  
    atom_object = read_cell(file_path, mag_element, pbc)
    atom_neigh = get_neighborlist(atom_object, rcut/2)
    mmx, mmy, mmz = get_mm(atom_object, magnetic_moment=mag_moment)
    dipole_z = calculate_MSA(atom_neigh, mmz)
    dipole_y = calculate_MSA(atom_neigh, mmy)
    dipole_x = calculate_MSA(atom_neigh, mmx)
    for i in range(len(dipole_x)):
        print('第{0}个原子:\n x方向: {1:10f} μeV, y方向: {2:.10f} μeV, z方向: {3:.10f} μeV'.format(i+1, dipole_x[i], dipole_y[i], dipole_z[i]))
        print(' out of plane - in plane: {0:.10f} μeV'.format(dipole_z[i]-min(dipole_x[i], dipole_y[i])))
    end_time = datetime.now()
    print(end_time-start_time)
```

运行代码，得到如下结果

![](https://pica.zhimg.com/80/v2-ecffde9fbc8dd2232b1208aa2c515a64_1440w.jpeg?source=d16d100b)

$E_{100}$= $E_x-E_z$ = −74.3μeV , E010=Ey−Ez=−64.2μeV ,与文献结果一致。\*参考文献Nanoscale, 2023, 15, 13402

文献中E100=-75μeV, E010=-64μeV![](https://pic1.zhimg.com/80/v2-6289aac74d886fb0954cc12375456f63_1440w.png?source=d16d100b)

脚本见MSA-calculate.py
