
# VASP 进阶性质计算

贡献者：张烜广

---

## 晶格热导率计算

参考教程：

1. [B 站视频](https://www.bilibili.com/video/BV1wo4y1U772?t=30.4) ；

2. [Phonopy 官网](https://phonopy.github.io/phonopy/index.html) ；

3. ShengBTE + thirdorder手册；

### 1. 高精度优化（VASP）

为保证优化收敛，可以精度逐渐调高，先同时优化晶格和原子位置用 `ISIF = 3`，后面只优化原子位置用 `ISIF = 2`。
逐渐增大 `EDIFFG` 的精度。INCAR参数设置如下：

```
EDIFFG = -1E-8
ISIF = 2
IBRION = 1
```

完美晶格优化具体步骤可以看 `wb/thermalconductivity/perfect` 文件夹下的 `opt2-opt5`（`opt1` 就是当前文件夹下的 VASP 计算文件）文件夹中的 `INCAR`。
（金刚石没那么难优化，其他的复杂体系可以按照这样的优化步骤。）

### 2. 有限位移法计算获得二阶力常数（Phonopy + VASP）

#### 2.1 对优化过后的晶胞扩胞

```shell
phonopy -d --dim='2 2 2'
```

#### 2.2 计算二阶力常数

对生成的 `POSCAR-\*` 文件全部进行静态计算，建立 `band.conf`。并依次执行一下命令：

```
phonopy -f vasprun.xml-00*
phonopy --dim='2 2 2' -c POSCAR-unitcell -p -s band.conf
phonopy-bandplot --gnuplot>222.dat
```

生成 `FORCE_CONSTANTS`。
本例中，完美晶格 `FORCE_CONSTANTS` 文件位置为 `wb/thermalconductivity/perfect/FC_2ND/FORCE_CONSTANTS`。

空位晶格同理。

计算声子谱，

1. 可以使用 `IBRION = 5/6/7/8`，`5` 和 `6` 是使用有限位移法进行计算，`7` 和 `8` 是使用紧束缚近似方法（DFPT）方法进行计算，紧束缚无法使用 `IVDW`。
2. 使用任何以上一个参数计算完成后，建立声子谱输入文件 `band.conf` 如下所示：

```
ATOM_NAME = Ge Se 
DIM = 5 5 1
BAND = 0.0 0.0 0.0  0.333 0.333 0.0  0.5 0.0 0.0  0.0 0.0 0.0
BAND_LABELS = G K M G
FORCE_CONSTANTS = READ
FC_SYMMETRY = .TRUE.
```

3. 使用命令 `phonopy -c POSCAR-unitcell -p -s band.conf` 直接输出声子谱pdf图片

声子谱详细计算信息参考 [网址](https://zhuanlan.zhihu.com/p/599295736) 。

> **ps**：
> band.conf中，设置 `FC_SYMMETRY=.TRUE.` 考虑了对称性，因此生成的 `FORCE_CONSTANT` 第一行两个数字（移动位置的原子数，扩胞原子数）并不相等。
> 将第二行命令更改为：
> ```
> phonopy --dim='2 2 2' -c POSCAR-unitcell band.conf --full-fc
> ```
> 生成不含对称性的 `FORCE_CONSTANT`，不然使用 ShengBTE 时会报错。

### 3. 三阶力常数矩阵计算（thirdorder + VASP）

#### 3.1 对优化过后的晶胞扩胞

```
thirdorder_vasp.py sow 2 2 2 -8
```

`2 2 2` 表示生成的超胞大小；`-8` 参数位为正数的时候表示截断半径（单位：nm）；`-8` 参数位为负数时表示第几近邻截断（此时应取负整数，例如-3表示第三近邻）。
对所有生成的 POSCAR 都进行高精度静态计算，INCAR 与有限位移法计算二阶力常数保持一致。

#### 3.2 计算三阶力常数

```shell
find job* -name vasprun.xml | sort -n | thirdorder_vasp.py reap 2 2 2 -8
```
生成`FORCE_CONSTANTS_3RD`。
本例中，完美晶格`FORCE_CONSTANTS_3RD`文件位置为 `wb/thermalconductivity/perfect/FC_3RD/FORCE_CONSTANTS_3RD`。

空位晶格同理。

> **ps**：
> `ChatGPT 4.0` 对于力常数的解释：
> 三阶力常数是指在固体物理学和材料科学中用来描述原子间相互作用随它们相对位置的三次变化而变化的参数。
> 这是一个高阶的概念，相对于更常见的二阶力常数而言，它涉及到更复杂的相互作用。
> 在晶格动力学中，二阶力常数涉及原子间作用力随两个原子之间的位移的二次变化，而三阶力常数则涉及作用力随位移的三次变化。具体来说：
> 
> 1.  二阶力常数：描述了当一个原子从其平衡位置微小位移时，作用于它的力如何随着这个位移线性变化（即力与位移成正比）。
> 
> 2.  三阶力常数：在此基础上，考虑原子位移更复杂的效应，即当原子从平衡位置发生位移时，作用于它的力不仅与位移成正比，还与位移的平方成比例。这意味着原子间的相互作用不仅取决于它们之间的距离，还取决于这些距离变化的方式。
> 
> 三阶力常数对于理解和描述非线性效应至关重要，比如非线性晶体振动、热导率的非线性行为以及材料的非线性声学性质。在热力学性质的计算、非线性声子相互作用的分析以及高压下物质行为的研究中，三阶力常数是一个重要的参数。通过这些高阶力常数，科学家能更准确地预测和解释复杂材料系统中的物理现象。

### 4. 晶格热导率计算（ShengBTE）

1. 将 `CONTROL`，`FORCE_CONSTANTS`，`FORCE_CONSTANTS_3RD` 放到一个文件夹下，将 `FORCE_CONSTANTS` 改为 `FORCE_CONSTANTS_2ND`。
2. `qsub ShengBTE.pbs` 提交任务。
3. 一共计算了 `100 K - 700 K`（间隔 `100 K`）七个温度下的晶格热导率。
	晶格热导率在每个温度文件夹中的 `TiK/`（`i` 代表温度）文件夹下的 `BTE.kappa_tensor` 文件最后一行。
	每列对应的方向顺序分别为 `XX`、`XY`、`XZ`、`YX`、`YY`、`YZ`、`ZX`、`ZY`、`ZZ`。
	取 `XX`、`YY`、`ZZ` 作为 `X`、`Y`、`Z` 方向的晶格热导率值。
	例如完美晶格在 300 K 下的晶格热导率的文件位置在 `wb/thermalconductivity/perfect/300K/T300K/BTE.kappa_tensor`。
	空位晶格同理。
4. 所有软件中的文件信息可以查阅软件说明书文件夹中的 ShengBTE-README.md。

---

## 电子热导率计算（BoltzTrap2 + VASP）

热导率 *$\kappa_{e}$* = 晶格热导率 *$\kappa_{e}$*(l) + 电子热导率 *$\kappa_{e}$*(e)

首先做高精度的结构优化和静态计算（使用的是和计算晶格热导率一样的 K 点分布）。
在静态计算文件夹中使用如下命令调用 `BoltzTrap2`:

```shell
btp2 -vv interpolate -m 3 wb/perfect/scf/case
```

生成 `interpolation.bt2` 文件。使用命令：

```
btp2 -vv integrate interpolation.bt2 200:700:50
```

得到的 `interpolation.trace` 文件记录了体系输运因子的信息。
`interpolation.trace` 中第 8 列 `kappae/tau0[W/(m*K*s)]` 是可以用来计算电子热导率的数据。

> **ps**
> 计算电子热导率的公式，以及公式中需要用到的参数：

---

## 热容计算（Phonopy）

参考教程：
1. [Phonopy 官网](https://phonopy.github.io/phonopy/index.html) 
2. [知乎教程](https://zhuanlan.zhihu.com/p/448370499)

首先，在计算二阶力常数的文件夹中建立 `mesh.conf` 文件。
使用命令 `phonopy -t mesh.conf -c POSCAR-unitcell` 计算热力学性质；
然后使用命令 `phonopy-propplot --gnuplot thermal_properties.yaml > thermal.dat` 生成热力学 `.dat` 文件方便作图。

本例中，完美晶格文件在 `wb/thermalconductivity/perfect/FC_2ND/thermal.dat`。
空位晶格同理。

---

## 热膨胀系数（Phonopy + VASP）

计算步骤：

1. 结构优化得到 CONTCAR 改成 POSCAR；

2. 改写缩放系数 0.95 - 1.05 建立 11 个 POSCAR；

3. 对 11 个 POSCAR 进行声子计算；

4. 计算热学性质并收集 `thermal_properties-{1..11}.yaml` 文件；

5. 提取 `e-v.dat` 文件；

6. `phonopy-qha` 计算热学性质。

完美晶格计算步骤见 `wb/thermalexpansivity/perfect/thermalexpansion.sh`；
热膨胀系数 `.dat` 文件位置为 `wb/thermalexpansivity/perfect/CP/thermal_expansion.dat`。
空位晶格同理。

---

## 杨氏模量，泊松比计算（Vaspkit + VASP）

教程：[公众号文章](https://mp.weixin.qq.com/s?__biz=MzI2OTQ4OTExOA==\&mid=2247487112\&idx=1\&sn=bfa1e8c7981b15e880cd8426c7f8a8ab\&chksm=eadec839dda9412f12c9020b189fd767a1d9d1bfbc1929822c9410bc8c24834a330b663b0038\&mpshare=1\&scene=1\&srcid=1104qk2jJQfIvkhhJR8u318F\&sharer_sharetime=1604472889425\&sharer_shareid=dbb0734284422c0aebeb989a7ec537c0#rd) ；

计算步骤：

1. 构建原胞；

2. 高精度优化；

3. 准备 `VPKIT.in` 文件，使用 `vaspkit \n200`生成发生应变的结构文件；
4. 使用步骤 2 中的 INCAR 对新生成的结构进行计算；
5. 将 `VPKIT.in` 文件第一行改为 `2`，使用 `vaspkit \n200` 计算；
6. 可以将弹性常数矩阵放入 [ELATE](https://progs.coudert.name/elate) 分析；

本例中，力学数据保存在：
完美晶格：`wb/Youngmodulus/perfect/opt/opt2/data.txt`
空位晶格：`wb/Youngmodulus/vacancy1/opt/data.txt`

---

## 光学性质计算（vaspkit+VASP）

教程：
1. [知乎帖子1](https://www.zhihu.com/question/287951798) ；
2. [知乎帖子2](https://zhuanlan.zhihu.com/p/669635613) ；
3. [公众号文章](https://mp.weixin.qq.com/s/or3KbvwAv_h3B68AJAB1fA) ；

计算步骤：

1. 结构优化；
2. 静态计算得到 WAVECAR 和 CHGCAR；
3. 使用静态计算的 WAVECAR 和 CHGCAR 进行光学计算；
   注意 INCAR 中 `NBANDS` 的取值一般为自洽计算 `OUTCAR` 的 `NBANDS` 的值的 2 - 3 倍；
4. 运行 `vaspkit \n711` 得到光学性质。

本例中，完美晶格结果在 `wb/optics/perfect/dielectricconstant`
光吸收系数：`ABSORPTION.dat`
能量损失谱：`ENERGY_LOSSSPECTRUM.dat`
消光系数：`EXTINCTION.dat`
复介电函数虚部：`IMAG.in`
复介电函数实部：`REAL.in`
反射系数：`REFLECTIVITY.dat`
折射系数：`REFRACTIVE.dat`

空位晶格同理，
输出文件 `ABSORB.dat`、`REFRACTIVE.dat`、`REFLECTIVITY.dat`、`EXTINCTION.dat` 和 `ENERGYLOSSSPECTRUM.dat`
依次为光吸收系数、折射系数、反射系数、消光系数和能量损失谱。

> **ps**：
> 1. 第一次计算后发现光学性质的曲线毛刺非常多，可能是光学中的 `NEDOS` 值没有取高，建立 `test2/` 文件夹设置 `NEDOS = 2000` 进行测试。
> 2. 精度仍然较差，建立 `test3/` 文件夹设置 `NEDOS = 10000` 进行测试。
> 3. 增大静态计算中的 K 点密度。生成新的 WAVECAR 和 CHGCAR。建立 `test4/` 文件夹进行测试。曲线变得比较平滑。所以光学性质参考 `test4/` 中的数据。

