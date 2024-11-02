# VASP 批量计算

本节介绍几个常用于 VASP 批量计算的脚本：`POSCAR.sh`、`link.sh`、`batchrunvaspkit.sh`、`batchrunvasp.pbs`、`check.sh` 和 `find.sh`，
并结合示例讲述如何使用脚本 `vasp2xyz.sh` 和 `vasp_to_xyz.py` 构建 `.xyz` 格式的势函数训练集。

## 生成结构

准备训练集之前，首先要明确做一个怎样的势函数？
通用势函数需要采样非常丰富的结构；而专注于研究某种性质的势函数，采样覆盖目标问题相关的势能面区域即可。

常规的势函数训练集包括对应温度区间的 AIMD 采样结构，结合一定的原子坐标微扰结构（微扰操作见下文）即可。

AIMD 的精度不用要求太高，所以如果 AIMD 的温度区间是低温、室温或中温时可以使用 VASP MLFF（VASP 版本大于等于 6.3）加速采样，具体操作见下文。
模拟过程中，AIMD 一般间隔 `0.1 ps` 对轨迹抽样，更小时间间隔内的结构之间具有强关联性，对增加训练集多样性无益，徒增训练集大小。

采样和微扰过后，需要重新进行高精度的单点能计算以获得高精度的原子受力和位力信息。
如果研究高压或辐照体系，采样之后的单点能计算可能需要非常高的精度（例如 `KSPACING = 0.1`、`SIGMA = 0.02`）来确保原子受力的准确。
此外，辐照相关的势函数还需要非常多的小团簇结构，以充分覆盖原子局域环境的剧烈变化。

### VASP MLFF 加速动力学采样

AIMD 采样时 INCAR 设置如下所示：
```plaintext
# 全局设置
NCORE = 12 # 设置合理的并行参数能加速模拟
ISTART = 0
ICHARG = 2
LWAVE = .F.
LCHARG = .F.
KSPACING = 2 # 使用较大的 KSPACING 加速模拟
KGAMMA = .T.
ISYM = 0 # AIMD 需要关闭对称性
# IVDW = 12 # 根据需要决定是否启用色散校正

# 自洽相关
ALGO = Normal # Fast 更快，但是 Normal 更稳定
NELM = 120
NELMIN = 4
EDIFF = 1E-05 # 使用较低的收敛标准加速模拟

# MD 相关，使用 NVT 系综
IBRION = 0
ISIF = 2
POTIM = 2.0 # 设置模拟步长
NSW = 5000 # 设置模拟步数
TEBEG = 300 # 设置模拟温度
TEEND = 300 # 设置模拟温度
MDALGO = 2
SMASS = 0
NBLOCK = 1

# MLFF 相关，不使用 MLFF 时全部注释即可
ML_MODE = train
ML_LMLFF = .T.
# ML_ISTART = 1 # 续算开关
# ML_CX = -0.1 # 恒温模拟时可以根据需要增大第一性原理计算比率
# ML_IALGO_LINREG = 3
# ML_RCUT2 = 6.0
# ML_RCUT1 = 6.0
# ML_CTIFOR = 1000
# ML_MB = 3000 # 模拟失败时根据LCONF适当增大
# ML_LBASIS_DISCARD = .TRUE.

# 精度相关
ENCUT = 400 # 减小 ENCUT 加速模拟
ISMEAR = 0
SIGMA = 0.05
PREC = N
LREAL = A
ADDGRID = .T.
```

设置 MLFF 时，常用命令如下：

```shell
grep STATUS ML_LOGFILE # 查看计算进程
tail -f ML_LOGFILE | grep STATUS # 实时监测进程
grep ERR ML_LOGFILE > mlff_err # 输出整个过程的误差变化
grep 'free  energy ML TOTEN' OUTCAR > mlff_ene # 得到整个过程的能量变化
```

如果想要续算 MLFF，在当前目录下执行
```shell
cp ML_ABN ML_AB
cp ML_FFN ML_FF
cp CONTCAR POSCAR
vi INCAR
```
并在 INCAR 中设置 `ML_ISTART = 1` 即可。

### 对平衡晶胞施加微扰

对平衡的晶胞施加微扰，一般包括 ±3% 的随机形变和每个原子 0.05、0.10、0.20 Å 的随机位移（如果模拟对象包含非晶相或液相，随机位移可以逐渐增大到1.0 Å）。
推荐使用更小的平衡晶胞（不超过64原子）来施加微扰，既不妨碍实现各种形变，又有效地控制了训练集中总的原子数目，是加速训练集构建行之有效的策略。
微扰后的结构可以使用 `screen_forces.py` 或类似脚本剔除不合理的结构（例如高嫩结构和原子受力特别大的结构）。

微扰脚本 `generate_random_perturbed_structures.py` 如下（仅供参考）：
```python
# Usage: python generate_random_perturbed_structures.py

import os
import numpy as np
from pymatgen.io.vasp.inputs import Poscar
from pymatgen.core.lattice import Lattice

def generate_random_perturbed_structures(filename, tolerance=[3, 3], n_struct=10):
    origin_struct = Poscar.from_file(filename, check_for_POTCAR=False).structure
    old_latt = origin_struct.lattice.as_dict(verbosity=1)
    old_latt = list(old_latt.values())[4:10]
    deform_struct = origin_struct
    delta, sigma = map(int, tolerance)
    Dlength = range(-delta, delta+1)
    Dangle = range(-sigma, sigma+1)
    lengths = old_latt[:3]
    angles = old_latt[3:]
    Pfile = 'perturb'
    if not os.path.exists(Pfile):
        os.makedirs(Pfile)
    for rnd in range(n_struct):
        lengths_loss = np.random.choice(Dlength, 3) / 100
        lengths_loss = np.array(lengths) * lengths_loss
        angle_loss = np.random.choice(Dangle, 3)
        loss = np.append(lengths_loss, angle_loss)
        Perturbed_lattice = np.array(old_latt) - loss
        a, b, c, alpha, beta, gamma = Perturbed_lattice
        new_latt = Lattice.from_parameters(a, b, c, alpha, beta, gamma)
        deform_struct.lattice = new_latt
        out_poscar = Poscar(deform_struct, comment=f'{rnd}_Perturbed structure')
        out_poscar.write_file(f'{Pfile}/perturbed_{rnd}.vasp')
    print(f"{n_struct} structures with random perturbation in lattice were generated! Bye!")

if __name__ == '__main__':
    file = 'POSCAR'  # 施加扰动的原始结构
    tolerance = [3, 3]  # 晶格矢量长度和晶格夹角扰动的容差百分比
    size = 50  # 生成扰动结构的数量
    generate_random_perturbed_structures(file, tolerance, size)
```

模拟对象无序度较高时，可以使用 RAG（Randomized Atomic-system Generator）或类似的方法迅速采样中高能区域势能面，但是这种方法可能引入非物理结构导致势能面失真，需要谨慎处理。
RAG 方法的具体细节参考 [Efficient Training of Machine Learning Potentials by a Randomized Atomic-System Generator](https://pubs.acs.org/doi/10.1021/acs.jpcb.0c05075) 一文。

---

## VASP 批量计算

假如手里已经有了筛选后的结构文件，例如 `iter0.xyz`，使用 Ovito 将其输出为一帧一帧的 `.vasp` 结构文件存放在文件夹 `iter0/` 下。
使用 `./POSCAR.sh iter0/` 将 `iter0/` 下的每一帧 `.vasp` 结构文件转为对应的子文件夹及子文件夹中的 `POSCAR` 文件。

`POSCAR.sh` 脚本如下：
```shell
# Usage: ./POSCAR.sh GeTe/
cd $1
for i in *
do
  name=$(basename $i .vasp)
  mkdir $name
  mv $i $name/POSCAR
done
cd ..
```

拷贝 `link.sh`、`batchrunvaspkit.sh`、`batchrunvasp.pbs`、`check.sh`、`find.sh` 和 `CONTCAR.sh` 到 `iter0/` 下，并在 `iter0/` 下准备批量计算用到的其它输入文件如 `INCAR`、`KPOINTS`、`POTCAR`。
这里给出一个单点计算的 `INCAR` 示例：
```plaintext
SYSTEM = GST

NCORE = 8
ISTART = 0
ICHARG = 2
LCHARG = .F.
LWAVE = .F.
KSPACING = 0.1
KGAMMA = .T.
ISYM = 0

ALGO = Fast
NELM = 120
NELMIN = 4
EDIFF = 1E-07

IBRION = -1
NSW = 0
EDIFFG = -1E-03

ENCUT = 600
ISMEAR = 0
SIGMA = 0.02
PREC = A
LREAL = A
ADDGRID = .T.
```

使用 `./link.sh INCAR KPOINTS POTCAR` 将输入文件链接到每一个子文件下。以 `INCAR` 为例，此时修改任何一个 `INCAR`，所有的 `INCAR` 都会改变。
`link.sh` 脚本如下：
```shell
# Usage: ./link.sh INCAR KPOINTS POTCAR
for i in "$@"
do
  for j in */
  do
    ln -s ../$i $j/$i
  done
done
```

如果结构文件的元素种类不同一，那么单点计算不能使用同样的 `POTCAR`；
类似的，如果结构文件的晶胞大小不一致，那么单点计算不能使用同样的 `KPOINTS`（使用 KSPACING 参数就没有这个烦恼）。
这是可以使用 `./batchrunvaspkit.sh` 批量调用 Vaspkit 生成需要的输入文件。

以生成 `KPOINTS` 为例，`batchrunvaspkit.sh` 脚本如下：
```shell
# Usage: ./batchrunvaspkit.sh
for i in */
do
  if [ -e $i/POSCAR ]
  then
    cd $i
    echo -e "102 \n2 \n0.04" | vaspkit
    cd $OLDPWD
  fi
done
```

输入文件齐全后，可以使用 `qsub batchrunvasp.pbs` 提交作业。
我个人的 `batchrunvasp.pbs` 脚本如下：
```shell
#PBS -N static
#PBS -l nodes=4:ppn=32
#PBS -l walltime=600:00:00
#PBS -q manycores
#PBS -V
#PBS -S /bin/bash

source /opt/intel/compilers_and_libraries_2018/linux/bin/compilervars.sh intel64
source /opt/intel/mkl/bin/mklvars.sh intel64
source /opt/intel/impi/2018.1.163/bin64/mpivars.sh

NP=`cat $PBS_NODEFILE | wc -l`
NN=`cat $PBS_NODEFILE | sort | uniq | tee /tmp/nodes.$$ | wc -l`

cat $PBS_NODEFILE > /tmp/nodefile.$$
cd $PBS_O_WORKDIR
ulimit -s unlimited

EXEC=/opt/software/vasp/vasp6.4.0/bin/vasp_std

for i in */
do
  if [ -e $i/vasprun.xml ]
  then
    echo 'skip ' $i >> ./tmp
  else
    cd $i
    mpirun -machinefile $PBS_NODEFILE -np $NP $EXEC > output
    time=`grep Elapsed OUTCAR`
    out=`tail -n 1 OUTCAR`
    echo $i >> ../tmp
    echo $time >> ../tmp
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
n=0
finished=0
for i in */
do
    ((n++))
    echo '-------------------' >> ./tmp
    if [ -e "${i}OSZICAR" ]; then
        ((finished++))
        step=$(grep "=" "${i}OSZICAR" | tail -1)
        conv=$(grep ":" "${i}OSZICAR" | tail -1 | awk '{print $2,$4}')
        gprra=$(grep "required" "${i}OUTCAR")
        echo "$i" "$step" >> ./tmp
        echo "$i" "$conv" >> ./tmp
        echo "$i" "$gprra" >> ./tmp
    else
        echo "$i" 'is still in line! ! !' >> ./tmp
    fi
done
echo "Total directories: $n"
echo "Unfinished directories: $((n - finished))"
```

考虑到每个人的磁盘空间是有限的，计算完成后删除所在目录下的 `CHG`、`CHGCAR`、`WAVECAR` 文件甚至 `POTCAR` 是很有必要的。
使用 `./find.sh CHG CHGCAR WAVECAR` 或者 `find . \( -name "CHG" -o -name "CHGCAR" -o -name "WAVECAR" \) -exec rm {} +` 删除不需要的文件。

```shell
# Usage: ./find.sh CHG CHGCAR WAVECAR
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

为了收集单点计算的数据，可以使用 `vasp_to_xyz.py` 脚本从 `vasprun.xml` 文件中提取 NEP 训练所需的 `.xyz` 格式的结构文件。
使用 `vasp2xyz.sh` 脚本可以批量提取并将结构文件合并到一个 `.xyz` 文件中。

`vasp_to_xyz.py` 脚本如下：

```python
# Usage: python vasp_to_xyz.py GST

from ase.io import read, write
import numpy as np
import sys
import os

def main(output_file='NEP-dataset.xyz', label='low', max_scf=120):
    os.system("find . -name vasprun.xml > xmllist")
    if os.path.exists('screen_tmp'):
        os.remove('screen_tmp')
    if os.path.exists(output_file):
        os.remove(output_file)

    with open('xmllist', 'r') as file_list:
        for line in file_list:
            xml = line.strip('\n')
            print(f"Processing file: {xml}")

            try:
                b = read(xml, index=":")
            except Exception:
                outcar_path = xml.replace("vasprun.xml", "OUTCAR")
                b = read(outcar_path, index=":")
                print(f"Fallback to: {outcar_path}")

            os.system(f"grep -B 1 E0 {xml.replace('vasprun.xml','OSZICAR')} | grep -E 'DAV|RMM' | awk '{{if($2>={max_scf}) print 0; else print 1}}' > screen_tmp")
            screen = np.loadtxt("screen_tmp")

            if screen.ndim == 0:
                screen = [screen]

            for ind, is_converged in enumerate(screen):
                if is_converged == 1:
                    xx, yy, zz, yz, xz, xy = -b[ind].calc.results['stress'] * b[ind].get_volume()
                    b[ind].info['virial'] = np.array([(xx, xy, xz), (xy, yy, yz), (xz, yz, zz)])
                    del b[ind].calc.results['stress']
                    b[ind].pbc = True
                    b[ind].info['config_type'] = label
                    write(output_file, b[ind], append=True)

    os.remove('screen_tmp')
    os.remove('xmllist')

if __name__ == "__main__":
    if len(sys.argv) == 2:
        label = sys.argv[1]
    else:
        label = 'low'
    main(label)
```

`vasp2xyz.sh` 脚本调用了上面的 `vasp_to_xyz.py` 脚本，具体代码如下：

```shell
# Usage: ./vasp2xyz.sh
for i in *
do
  if [ -e $i/vasprun.xml ]
  then
    cd $i
    python vasp_to_xyz.py GST # 可以指定标签
    echo $i
    cat NEP-dataset.xyz >> ../NEP-dataset.xyz
    cd $OLDPWD
  fi
done
```
