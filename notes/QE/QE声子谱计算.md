
# QE 声子谱计算

## 结构优化

声子谱计算需要非常高的计算精度，并确保初始结构已充分优化以减少虚频的可能性。

```
&CONTROL
 calculation = 'vc-relax'
 verbosity = 'high'
 nstep = 100
 tstress = .true.
 tprnfor = .true.
 outdir = './outdir'
 prefix = 'Si'
 etot_conv_thr = 1.0d-7
 forc_conv_thr = 1.0d-5
 disk_io = 'low'
 pseudo_dir = '/home/changruiwang-ICME/UPF/SSSP/'
/
&SYSTEM
 ibrav = 2
 celldm(1) = 10.410909236
 nat = 2
 ntyp = 1
 nbnd = 8
 ecutwfc = 50.0
 ecutrho = 200.0
 occupations = 'fixed'
 smearing = 'gaussian'
 degauss = 0.01
/
&ELECTRONS
 electron_maxstep = 300
 scf_must_converge = .false.
 conv_thr = 1.0d-9
 mixing_mode = 'plain'
 mixing_beta = 0.8d0
 diagonalization = 'david'
/
&IONS
 ion_dynamics = 'bfgs'
 trust_radius_min = 1.0d-3
 trust_radius_ini = 0.5d0
/
&CELL
 cell_dynamics = 'bfgs'
 press = 0.0
 press_conv_thr = 0.1
 cell_dofree = 'all'
/
ATOMIC_SPECIES
 Si  28.0855  Si.pbe-n-rrkjus_psl.1.0.0.UPF
ATOMIC_POSITIONS (alat)
 Si 0.00 0.00 0.00
 Si 0.25 0.25 0.25
K_POINTS (automatic)
  12 12 12 0 0 0
```

运行 `pw.x < 1_vc_relax.inp > 1_vc_relax.out` 即可。

计算结束后，运行

```
awk '/Begin final coordinates/,/End final coordinates/{print $0}' output
```

查看 vc-relax 的结果。

---

## 自洽计算

vc-relax 计算之后需要再做一次自洽计算，否则后续计算会出错。
修改 `calculation = 'vc-relax'` 为 `calculation = 'scf'`，然后运行 `pw.x < 2_scf.inp > 2_scf.out` 即可。
当 `calculation = 'scf'` 时，`pw.x` 会自动忽略 `IONS` 和 `CELL` 等不相干的模块。

---

## 声子计算

自洽计算完成后，需要在倒空间进行动力学矩阵的计算。这一步耗时最长。

### 单个 q 点的计算

```
phonons of Si at Gamma
 &INPUTPH
  amass = 28.0855 # 每种原子类型的原子质量 [amu]
  outdir = './outdir' # 包含输入、输出和临时文件的目录
  prefix = 'Si' # 输入/输出文件的前缀
  tr2_ph = 1.0d-14 # 自洽收敛的阈值
  fildyn = 'Si.dynG' # 用于保存动力学矩阵的文件。
  epsil = .true. # 计算系统的宏观介电常数，不适用于金属体系
 /
0.0 0.0 0.0
```

运行 `ph.x < 3_ph.inp > 3_ph.out` 即可。

### 声子色散计算

```
phonons of Si q-grid
 &INPUTPH
  amass = 28.0855
  outdir = './outdir'
  prefix = 'Si'
  tr2_ph = 1.0d-12
  fildyn = 'Si444.dyn'
  ldisp = .true. # 控制计算声子色散
  nq1 = 4 # Monkhorst-Pack 网格参数
  nq2 = 4
  nq3 = 4
 /
```

运行 `ph.x < 3_ph.inp > 3_ph.out` 即可。

---

## 后处理

`ph.x` 计算结束后，使用 `q2r.x` 读取由 ph.x 代码生成的 q 点网格的力常数矩阵 C(q)，做傅里叶变换，并在实空间计算相应的原子间力常数 (IFC) C(R)。

输入 `q2r.inp` 如下，再运行 `q2r.x < 4_q2r.inp > 4_q2r.out` 即可。

```
&INPUT
   fildyn = 'Si444.dyn' # 输入的动力学矩阵文件
   flfrc = 'Si444.fc' # 输出的原子间力常数文件
   zasr = 'crystal' # 'simple' 适合简单晶体系统或粗略计算，'crystal' 适合大多数晶体计算
 /
```

在输出文件 `4_q2r.out` 中需要有以下内容：
```
fft-check success (sum of imaginary terms < 10^-12)
```

然后使用 `matdyn.x` 在超胞基础上计算任意 q 点的声子频率。
将以下内容保存为 `5_matdyn.inp`，运行 `matdyn.x < 5_matdyn.inp > 5_matdyn.out`。

```
&INPUT
  flfrc = 'Si444.fc' # 由 q2r.x 生成的力常数文件
  asr = 'crystal' # Acoustic Sum Rule 类型
  flfrq = 'Si.freq' # 输出的声子频率文件
  q_in_band_form = .true. # q 点以带的形式给出，指定起点和终点即可。
  !q_in_cryst_coord=.true. # q 点以坐标的形式给出
/
5
  gG 30
  X 30
  W 30
  K 30
  gG 30
```

最后，使用 `plotband.x` 来处理输出的声子频率文件。
将以下内容保存为 `6_plotband.inp` 文件，然后运行 `plotband.x < 6_plotband.inp > 6_plotband.out`，生成的 `freq.plot` 文件可以直接导入 Origin 做图。

```
Si.freq # 读入的数据文件名
0 500 # 能量范围下限和上限
freq.plot # 输出声子色散带数据
freq.ps # 输出可视化的图片
0.0 # Fermi 能，能带做图相关，这里不影响
50.0 0.0 # 能量度量间隔和能量平移大小
```
