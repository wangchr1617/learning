# CP2K 批量计算

本节介绍几个常用于 CP2K 批量计算的脚本：`CIF.sh`、`link.sh`、`batchruncp2k.pbs`、`check.sh` 和 `find.sh`，
并结合示例讲述如何使用脚本 `cp2k2dp.py` 和 `dp2xyz.py` 构建 `.xyz` 格式的势函数训练集。

## 生成结构

准备训练集之前，首先要明确做一个怎样的势函数？
通用势函数需要采样非常丰富的结构；而专注于研究某种性质的势函数，采样覆盖目标问题相关的势能面区域即可。

常规的势函数训练集包括对应温度区间的 AIMD 采样结构，结合一定的原子坐标微扰结构（微扰操作见下文）即可。

AIMD 的精度不用要求太高，例如可以使用 TZVP 基组、100 Ry 截断、单 Gamma 点加速采样。
使用 CP2K 进行 AIMD 模拟时，要注意开启 MOTION/PRINT/CELL 和 MOTION/PRINT/FORCES，如下所示：

```plaintext
&MOTION
  &MD
    ……
  &END MD
  &PRINT
    &TRAJECTORY
      &EACH
        MD     1
      &END EACH
    &END TRAJECTORY
    &VELOCITIES
      &EACH
        MD     0
      &END EACH
    &END VELOCITIES
    &FORCES
      &EACH
        MD     1
      &END EACH
    &END FORCES
    &CELL
      &EACH
        MD     1
      &END EACH
    &END CELL
    &RESTART
      BACKUP_COPIES 0
      &EACH
        MD  1
      &END EACH
    &END RESTART
    &RESTART_HISTORY
      &EACH
        MD  0 # 另一种策略是指定输出 restart 文件的频率，然后使用 cp2k2vasp.sh 将 restart 文件转成 POSCAR
      &END EACH
    &END RESTART_HISTORY
  &END PRINT
&END MOTION
```

否则无法使用 `cp2k2xyz.py` 脚本提取 NEP 训练所需的 `.xyz` 格式的结构文件。
`cp2k2xyz.py` 脚本如下：

```python
# Usage: python cp2k2xyz.py ./

import sys
import os

def get_cp2k_filenames(directory):
    position_file = None
    force_file = None
    cell_file = None
    for filename in os.listdir(directory):
        if filename.endswith(".xyz"):
            if "-pos-" in filename:
                position_file = filename
            elif "-frc-" in filename:
                force_file = filename
        elif filename.endswith(".cell"):
            cell_file = filename
    if position_file is None or force_file is None or cell_file is None:
        print(f"请检查目录 {directory} 中的文件。")
        raise FileNotFoundError("找不到必要的文件: 位置文件, 力文件或晶胞文件缺失。")
    return position_file, force_file, cell_file

def convert_cp2k_to_xyz(directory, filenames=None, output_filename=None):
    if filenames is None:
        position_file, force_file, cell_file = get_cp2k_filenames(directory)
    else:
        position_file, force_file, cell_file = filenames
    os.makedirs("CP2K2XYZ", exist_ok=True)
    if output_filename is None:
        output_filename = "merged_cp2k.xyz"
    position_path = os.path.join(directory, position_file)
    force_path = os.path.join(directory, force_file)
    cell_path = os.path.join(directory, cell_file)
    output_path = os.path.join("CP2K2XYZ", output_filename)
    print("输入文件:", position_path, force_path, cell_path)
    print("输出文件:", output_path)
    with open(position_path, 'r') as pf, open(force_path, 'r') as ff, \
         open(cell_path, 'r') as cf, open(output_path, 'w') as of:
        cf.readline()
        while True:
            position_header = pf.readline()
            force_header = ff.readline()
            if not position_header or not force_header:
                break
            num_atoms = int(position_header.strip().split()[0])
            of.write(f"{num_atoms}\n")
            force_info_line = ff.readline()
            pf.readline()
            energy = float(force_info_line.strip().split("E =")[-1]) * 27.211386245988
            cell_line = cf.readline().strip().split()
            lattice = " ".join(cell_line[2:11])
            of.write(f"energy={energy:.10f} config_type=cp2k2xyz pbc=\"T T T\" ")
            of.write(f"Lattice=\"{lattice}\" Properties=species:S:1:pos:R:3:forces:R:3\n")
            for _ in range(num_atoms):
                position_line = pf.readline().strip().split()
                force_line = ff.readline().strip().split()
                if len(position_line) < 4 or len(force_line) < 4:
                    break
                force_x = float(force_line[1]) * 51.42206747632590000
                force_y = float(force_line[2]) * 51.42206747632590000
                force_z = float(force_line[3]) * 51.42206747632590000
                of.write(f"{position_line[0]} {position_line[1]} {position_line[2]} {position_line[3]} ")
                of.write(f"{force_x:.10f} {force_y:.10f} {force_z:.10f}\n")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("用法: python cp2k2xyz.py <目录路径>")
        sys.exit(1)

    directory_path = sys.argv[1]
    convert_cp2k_to_xyz(directory_path)
```

时间步长设置为 2 fs，间隔 `0.1 ps` 对轨迹抽样即可。
更小时间间隔内的结构之间具有强关联性，对增加训练集多样性无益，徒增训练集大小。

采样和微扰过后，需要重新进行高精度的单点能计算以获得高精度的原子受力。
值得注意的是，CP2K 单点计算得到的位力信息不是很准，所有有时引入位力信息对势函数训练反而有害无益。

### 对平衡晶胞施加微扰

对平衡的晶胞施加微扰，一般包括 ±3% 的随机形变和每个原子 0.05、0.10、0.20 Å 的随机位移（如果模拟对象包含非晶相或液相，随机位移可以逐渐增大到1.0 Å）。
推荐使用更小的平衡晶胞（不超过64原子）来施加微扰，既不妨碍实现各种形变，又有效地控制了训练集中总的原子数目，是加速训练集构建行之有效的策略。
微扰后的结构可以使用 `screen_forces.py` 或类似脚本剔除不合理的结构（例如高嫩结构和原子受力特别大的结构）。

模拟对象无序度较高时，可以使用 RAG（Randomized Atomic-system Generator）或类似的方法迅速采样中高能区域势能面，但是这种方法可能引入非物理结构导致势能面失真，需要谨慎处理。
RAG 方法的具体细节参考 [Efficient Training of Machine Learning Potentials by a Randomized Atomic-System Generator](https://pubs.acs.org/doi/10.1021/acs.jpcb.0c05075) 一文。

---

## CP2K 批量计算

假如手里已经有了筛选后的结构文件，例如 `iter0.xyz`，使用 Ovito 将其输出为一帧一帧的 `.cif` 结构文件存放在文件夹 `iter0/` 下。
使用 `./CIF.sh iter0/` 将 `iter0/` 下的每一帧 `.cif` 结构文件转为对应的子文件夹及子文件夹中的 `input.cif` 文件。

`CIF.sh` 脚本如下：

```shell
# Usage: ./CIF.sh GeTe/
cd $1
for i in *
do
  name=$(basename $i .cif)
  mkdir $name
  mv $i $name/input.cif
done
cd ..
```

拷贝 `link.sh`、`batchruncp2k.pbs`、`check.sh` 和 `find.sh` 到 `iter0/` 下，并在 `iter0/` 下准备批量计算用到的输入文件 `cp2k.inp`。
这里给出一个单点计算的 `cp2k.inp` 示例：

```plaintext
#Generated by Multiwfn
&GLOBAL
  PROJECT cp2k
  PRINT_LEVEL MEDIUM # 一定要设置成 MEDIUM
  RUN_TYPE ENERGY_FORCE
  PREFERRED_DIAG_LIBRARY SL # 用于解决 ELPA 库缺失的问题
&END GLOBAL

&FORCE_EVAL
  METHOD Quickstep
  &SUBSYS
    &TOPOLOGY # 使用读取文件的方式来获取原子坐标
      COORD_FILE_FORMAT CIF 
      COORD_FILE_NAME input.cif
    &END TOPOLOGY
    &CELL # 使用读取文件的方式来获取晶格参数
      CELL_FILE_FORMAT CIF
      CELL_FILE_NAME input.cif
      PERIODIC XYZ
    &END CELL
    &KIND Ge
      ELEMENT Ge
      BASIS_SET TZVP-MOLOPT-PBE-GTH-q4
      POTENTIAL GTH-PBE
    &END KIND
    &KIND Te
      ELEMENT Te
      BASIS_SET TZVP-MOLOPT-PBE-GTH-q6
      POTENTIAL GTH-PBE
    &END KIND
  &END SUBSYS

  &DFT
    BASIS_SET_FILE_NAME  BASIS_MOLOPT_UZH
    POTENTIAL_FILE_NAME  POTENTIAL
    CHARGE    0 
    MULTIPLICITY    1 
    &KPOINTS # 根据原子数合理地设置 K 点
      SCHEME MONKHORST-PACK  2  2  2
    &END KPOINTS
    &QS
      EPS_DEFAULT 1.0E-14
    &END QS
    &POISSON
      PERIODIC XYZ
      PSOLVER PERIODIC
    &END POISSON
    &XC
      &XC_FUNCTIONAL PBE
      &END XC_FUNCTIONAL
    &END XC
    &MGRID
      CUTOFF  400
      REL_CUTOFF  55
    &END MGRID
    &SCF
      MAX_SCF 128
      EPS_SCF 1.0E-07
      &DIAGONALIZATION
        ALGORITHM STANDARD
      &END DIAGONALIZATION
      &MIXING
        METHOD BROYDEN_MIXING 
        ALPHA 0.4
        NBROYDEN 8 
      &END MIXING
      &SMEAR
        METHOD FERMI_DIRAC
        ELECTRONIC_TEMPERATURE 300 
      &END SMEAR
      ADDED_MOS 96 # 根据原子数合理地设置虚拟轨道数，一般设置成原子总数的 1/2 即可
      &PRINT
        &RESTART OFF
        &END RESTART
      &END PRINT
    &END SCF
  &END DFT
  &PRINT
    &FORCES ON # 输出原子受力
    &END FORCES
  &END PRINT
&END FORCE_EVAL
```

使用 `./link.sh cp2k.inp` 将输入文件链接到每一个子文件下。此时修改任何一个 `cp2k.inp`，所有的 `cp2k.inp` 都会改变。
`link.sh` 脚本如下：

```shell
# Usage: ./link.sh cp2k.inp
for i in "$@"
do
  for j in */
  do
    ln -s ../$i $j/$i
  done
done
```

如果结构文件的元素种类不同一，那么单点计算不能使用同样的 `cp2k.inp`。
建议根据晶胞原子数将相同原子数的结构放到同一个文件夹下使用相同的 `cp2k.inp` 进行单点计算。

输入文件齐全后，可以使用 `qsub batchruncp2k.pbs` 提交作业。
我个人的 `batchruncp2k.pbs` 脚本如下：

```shell
#PBS -N static
#PBS -l nodes=3:ppn=24
#PBS -l walltime=600:00:00
#PBS -q manycores
#PBS -V
#PBS -S /bin/bash

NP=`cat $PBS_NODEFILE | wc -l`
NN=`cat $PBS_NODEFILE | sort | uniq | tee /tmp/nodes.$$ | wc -l`

cat $PBS_NODEFILE > /tmp/nodefile.$$
cd $PBS_O_WORKDIR
ulimit -s unlimited

module load apptainer/1.0.0
EXEC="apptainer exec /home/jingjinghu-ICME/Software/cp2k-2024/cp2k-2024.sif cp2k.psmp"
export OMP_NUM_THREADS=1

for i in */
do
  if [ -e $i/cp2k.out ]
  then
    echo 'skip ' $i >> ./tmp
  else
    cd $i
    mpirun -machinefile $PBS_NODEFILE -np $NP $EXEC -i cp2k.inp 1>cp2k.out 2>cp2k.err
    out=$(grep "SCF run converged in" cp2k.out)
    echo $i >> ../tmp
    echo $out >> ../tmp
    cd ..
  fi
done
```

计算过程中，每一个子任务完成都会在 `./tmp` 文件中输出对应信息。计算完成后，使用 `./check.sh` 检查收敛情况。
`check.sh` 脚本如下：

```shell
# Usage: ./check.sh
rm -rf ./tmp

ntot=0
nfinished=0

for i in */
do
    ((ntot++))
    echo '-------------------' >> ./tmp
    if [ -e "${i}cp2k.out" ]; then
        ((nfinished++))
        conv=$(grep "SCF run converged in" "${i}cp2k.out")
        warn=$(grep "The number of warnings" "${i}cp2k.out")
        echo "$i" "$conv" >> ./tmp
        echo "$warn" >> ./tmp
    else
        echo "$i" '! ! !' >> ./tmp
    fi
done

echo "Total directories: $ntot"
echo "Unfinished directories: $((ntot - nfinished))"
```

考虑到每个人的磁盘空间是有限的，计算完成后删除所在目录下的 `.wfn` 文件（如果有设置输出的话）是很有必要的。
使用 `./find.sh *.wfn` 删除不需要的文件。
`find.sh` 脚本如下：

```shell
# Usage: ./find.sh *.wfn
for i in "$@"
do
  find ./ -name "$i" >> name.sh
done
sed -i 's/^/rm /' name.sh
chmod u+x name.sh
./name.sh
rm name.sh
```

---

## 单点数据收集

为了收集单点计算的数据，可以使用 `cp2k2dp.py` 脚本从 `cp2k.out` 文件中提取 DeePMD 训练所需的结构文件。
然后，使用 `dp2xyz.py` 脚本将 DeePMD 格式的结构文件批量转化成 NEP 训练所需的 `.xyz` 文件。

`cp2k2dp.py` 脚本如下：

```python
# Usage: python cp2k2dp.py

import dpdata
import numpy as np

data = dpdata.MultiSystems.from_dir('./', file_name='cp2k.out', fmt='cp2k/output', type_map=['Ge','Te'])
print(data)
print(data.systems)
data.to_deepmd_npy('../01_fp/',fmt='deepmd-npy')
```
注意指定路径和 `type_map`。

`dp2xyz.py` 脚本如下：

```python
# Usage: python dp2xyz.py 01_fp/
"""
    Purpose:
        Convert deepmd input file format to xyz.
    Ref:
        dpdata: https://github.com/deepmodeling/dpdata
    Run:
        python deep2xyz.py deepmd
"""

import os
import sys
import glob
import numpy as np

ELEMENTS=['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', \
         'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag',\
         'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb',\
         'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U', 'Np', \
         'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr']

def vec2volume(cells):
    va = cells[:3]
    vb = cells[3:6]
    vc = cells[6:]
    return np.dot(va, np.cross(vb,vc))

def cond_load_data(fname) :
    tmp = None
    if os.path.isfile(fname) :
        tmp = np.load(fname)
    return tmp

def load_type(folder, type_map=None) :
    data = {}
    data['atom_types'] = np.loadtxt(os.path.join(folder, 'type.raw'), ndmin=1).astype(int)
    ntypes = np.max(data['atom_types']) + 1
    data['atom_numbs'] = []
    for ii in range (ntypes) :
        data['atom_numbs'].append(np.count_nonzero(data['atom_types'] == ii))
    data['atom_names'] = []
    if os.path.isfile(os.path.join(folder, 'type_map.raw')) :
        with open(os.path.join(folder, 'type_map.raw')) as fp:
            my_type_map = fp.read().split()
    elif type_map is not None:
        my_type_map = type_map
    else:
        my_type_map = []
        for ii in range(ntypes) :
            my_type_map.append('Type_%d' % ii)
    assert(len(my_type_map) >= len(data['atom_numbs']))
    for ii in range(len(data['atom_numbs'])) :
        data['atom_names'].append(my_type_map[ii])

    return data

def load_set(folder) :
    cells = np.load(os.path.join(folder, 'box.npy'))
    coords = np.load(os.path.join(folder, 'coord.npy'))
    eners  = cond_load_data(os.path.join(folder, 'energy.npy'))
    forces = cond_load_data(os.path.join(folder, 'force.npy'))
    virs   = cond_load_data(os.path.join(folder, 'virial.npy'))
    return cells, coords, eners, forces, virs

def to_system_data(folder):
    data = load_type(folder)
    data['orig'] = np.zeros([3])
    data['docname'] = folder
    sets = sorted(glob.glob(os.path.join(folder, 'set.*')))
    all_cells = []
    all_coords = []
    all_eners = []
    all_forces = []
    all_virs = []
    for ii in sets:
        cells, coords, eners, forces, virs = load_set(ii)
        nframes = np.reshape(cells, [-1,3,3]).shape[0]
        all_cells.append(np.reshape(cells, [nframes,3,3]))
        all_coords.append(np.reshape(coords, [nframes,-1,3]))
        if eners is not None:
            eners = np.reshape(eners, [nframes])
        if eners is not None and eners.size > 0:
            all_eners.append(np.reshape(eners, [nframes]))
        if forces is not None and forces.size > 0:
            all_forces.append(np.reshape(forces, [nframes,-1,3]))
        if virs is not None and virs.size > 0:
            all_virs.append(np.reshape(virs, [nframes,9]))
    data['frames'] = nframes
    data['cells'] = np.concatenate(all_cells, axis = 0)
    data['coords'] = np.concatenate(all_coords, axis = 0)
    if len(all_eners) > 0 :
        data['energies'] = np.concatenate(all_eners, axis = 0)
    if len(all_forces) > 0 :
        data['forces'] = np.concatenate(all_forces, axis = 0)
    if len(all_virs) > 0:
        data['virials'] = np.concatenate(all_virs, axis = 0)
    if os.path.isfile(os.path.join(folder, "nopbc")):
        data['nopbc'] = True
    return data

def read_multi_deepmd(folder):
    data_multi = {}
    list_dir = []
    for dirpath, filedir, filename in os.walk(folder):
        if 'type.raw' in filename:
            list_dir.append(dirpath)
    for i, fi in enumerate(list_dir):
        ifold = fi
        idata = to_system_data(ifold)
        if 'virials' in idata and len(idata['virials']) == idata['frames']:
            idata['has_virial'] = np.ones((idata['frames']), dtype=bool)
        else:
            idata['has_virial'] = np.zeros((idata['frames']), dtype=bool)
        data_multi[i] = idata
    nframes = np.sum([data_multi[i]['frames'] for i in data_multi])
    data = {}
    data['nframe'] = nframes
    data['docname'] = ['' for i in range(nframes)]
    data['atom_numbs'] = np.zeros((data['nframe']))
    data['has_virial'] = np.zeros((data['nframe']))
    data['energies'] = np.zeros((data['nframe']))
    data['virials'] = np.zeros((data['nframe'], 9))
    data['cells'] = np.zeros((data['nframe'], 9))
    data['volume'] = np.zeros((data['nframe']))
    data['atom_names'] = {}
    data['atom_types'] = {}
    data['coords'] = {}
    data['forces'] = {}
    ifr = -1
    for i in data_multi:
        atom_names = [data_multi[i]['atom_names'][j] for j in data_multi[i]['atom_types']]
        for j in range(data_multi[i]['frames']):
            ifr += 1
            data['atom_numbs'][ifr] = len(data_multi[i]['atom_types'])
            data['has_virial'][ifr] = data_multi[i]['has_virial'][j]
            data['energies'][ifr] = data_multi[i]['energies'][j]
            if data['has_virial'][ifr]:
                data['virials'][ifr] = data_multi[i]['virials'][j]
            data['cells'][ifr] = np.reshape(data_multi[i]['cells'][j],9)
            data['volume'][ifr] = vec2volume(data['cells'][ifr])
            data['atom_names'][ifr] = atom_names
            data['atom_types'][ifr] = data_multi[i]['atom_types']
            data['coords'][ifr] = data_multi[i]['coords'][j]
            data['forces'][ifr] = data_multi[i]['forces'][j]
            data['docname'][ifr] = data_multi[i]['docname']
    return data

def check_data(data):
    print('Nframes:', data['nframe'])
    for i in range(data['nframe']):
        print(i, data['docname'][i])
        print('    atom_numbs', int(data['atom_numbs'][i]), end=' ')
        print('atom_types', len(data['atom_types'][i]))
        print('    coords', len(data['coords'][i]), end=' ')
        print('forces', len(data['forces'][i]))

def dump_nep(folder, data, nep_version=3):
    os.makedirs(folder, exist_ok=True)
    fout = open(os.path.join(folder, 'train.in'), 'w')
    outstr = str(data['nframe']) + '\n'
    for i in range(data['nframe']):
        outstr=outstr+str(int(data['atom_numbs'][i]))+' '+str(int(data['has_virial'][i]))+'\n'
    fout.write(outstr)
    for i in range(data['nframe']):
        outstr = ''
        if data['has_virial'][i]:
            outstr=outstr+str(data['energies'][i])+' '
            outstr=outstr+' '.join(map(str, data['virials'][i][[0, 4, 8, 1, 5, 2]]))+'\n'
        else:
            outstr=outstr+str(data['energies'][i])+'\n'
        outstr=outstr+' '.join(map(str, data['cells'][i]))+'\n'
        for j in range(int(data['atom_numbs'][i])):
            if nep_version == 1:
                ijname=data['atom_names'][i][j]
                ijanum=ELEMENTS.index(data['atom_names'][i][j]) + 1
                outstr=outstr+str(int(ijanum))+' '
            elif nep_version == 2:
                ijtype=data['atom_types'][i][j]
                outstr=outstr+str(int(ijtype))+' '
            elif nep_version == 3:
                ijname=data['atom_names'][i][j]
                outstr=outstr+ijname+' '
            else:
                raise "Errors with wrong <nep_version> para."
            outstr=outstr+' '.join(map(str, data['coords'][i][j]))+' '
            outstr=outstr+' '.join(map(str, data['forces'][i][j]))+'\n'
        fout.write(outstr)
    fout.close()

def dump_xyz(folder, data, outxyz="NEP-dataset.xyz"):
    os.makedirs(folder, exist_ok = True)
    Out_string = ""
    for i in range(data['nframe']):
        Out_string += str(int(data['atom_numbs'][i])) + "\n"
        myvirial = data['virials'][i]
        Out_string += "energy=" + str(data['energies'][i]) + " "
        Out_string += "config_type=dp2xyz "
        Out_string += "pbc=\"T T T\" "
        if data['has_virial'][i]:
            Out_string += "virial=\"" + " ".join(list(map(str, myvirial))) + "\" "
        Out_string += "Lattice=\"" + " ".join(list(map(str, data['cells'][i]))) + "\" "
        Out_string += f"volume={data['volume'][i]} "
        Out_string += "Properties=species:S:1:pos:R:3:forces:R:3\n"
        for j in range(int(data['atom_numbs'][i])):
            Out_string += data['atom_names'][i][j] + " "
            Out_string += " ".join(list(map(str, data['coords'][i][j]))) + " "
            Out_string += " ".join(list(map(str, data['forces'][i][j]))) + "\n"
    fo = open(os.path.join(folder, outxyz), 'w')
    fo.write(Out_string)
    fo.close()

def main():
    instr = sys.argv[1]
    data = read_multi_deepmd('./'+instr)
    dump_xyz('./XYZ', data, outxyz="NEP-dataset.xyz")

if __name__ == "__main__":
    main()
```
