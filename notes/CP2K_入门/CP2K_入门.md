# CP2K_入门

CP2K 的输入文件一般命名为 `cp2k.inp`，其中关键词以 `section` 和 `subsection` 的形式，一环套一环，如下图所示。

<div align="left">
<img src="./figures/fig_001.png" width = "50%" />
</div>

每一个 `section` 都是以 `&section` 开头，以 `&END section` 结尾，顺序随意，但是嵌套不能乱。
在 `section` 中，每行只写一个关键词，关键词后面接参数，CP2K 对大小写和空格不敏感。

常用的 `section` 包括 `GLOBAL`、`FORCE_EVAL`、`MOTION` 等，下面将逐个介绍它们。

---

## GLOBAL

`GLOBAL` 控制 CP2K 的全局参数，`GLOBAL` 的设置如下所示：
```
&GLOBAL
  PROJECT cp2k
  RUN_TYPE GEO_OPT
  PRINT_LEVEL LOW
&END GLOBAL 
```

`PROJECT` 指定计算任务名，个人习惯是统一设置成 `cp2k`，方便脚本后处理。

`RUN_TYPE` 指定 CP2K 计算任务的类型，包括
1. `MD`，分子动力学模拟；
2. `ENERGY`，单点能计算；
3. `ENERGY_FORCE`，计算能量和原子受力；
4. `GEO_OPT`，几何优化；
5. `CELL_OPT`，晶胞优化，常与 `GEO_OPT` 联用；
6. `BAND`，NEB计算；

`PRINT_LEVEL` 控制输出详略，一般设置成 `LOW` 或 `MEDIUM` 即可。

CP2K 中计算能量需要设置 `GLOBAL/RUN_TYPE` 为 `ENERGY`；
如果还要计算受力，需要设置 `GLOBAL/RUN_TYPE` 为 `ENERGY_FORCE`，
并设置 `FORCE_EVAL/PRINT/FORCES` 为 `ON`。
如果要将受力信息存储到文件中，则需要设定 `FROCE_EVAL/PRINT/FORCES/FILENAME` 文件名；
否则受力信息将会打印到 `*.out` 文件中。

注意，在 CP2K 中，默认的能量和距离单位分别是 `Hartree` (Hartree 原子单位) 和 `Bohr` (波尔半径)。
如果需要将它们转换成常用的 `eV`（电子伏特）和 `Å`（埃）单位，可以使用以下换算关系：
```
1 Hartree = 27.2114 eV
1 Bohr = 0.529177 Å
```
反之：
```
1 eV = 0.0367493 Hartree
1 Å = 1.88973 Bohr
```

---

## FORCE_EVAL

`FORCE_EVAL` 控制能量和原子受力（类似于 VASP 中的电子步），是 CP2K 中最重要的 `section`。
以 `SiSb` 的计算为例，一个典型的 `FORCE_EVAL`（使用 `DIAG` 算法）如下所示：
```
&FORCE_EVAL
  METHOD QUICKSTEP # FIST 是经典 MD，QUICKSTEP 是 AIMD；
  # STRESS_TENSOR ANALYTICAL # 对于晶胞体积/形状变化的计算需要打开这个参数；  
  &SUBSYS
    &CELL
	  # 下面两种指定晶胞参数的方法是等价的；
      A 18.62159700 0.000000000 0.000000000
      B 0.000000000 18.62159700 0.000000000
      C 0.000000000 0.000000000 18.62159700
      # ABC 18.62159700 18.62159700 18.62159700
      # ALPHA_BETA_GAMMA 90 90 90

      PERIODIC XYZ # 默认 XYZ，可选 X、Y、Z、XY、YZ、XZ、XYZ 和 NONE；
    &END CELL
    # &COORD # 指定原子坐标，更推荐使用下面的 TOPOLOGY 方法；
      # Si 0.000000000 0.000000000 0.000000000
      # Si 0.000000000 2.715348700 2.715348700
      # …
      # Sb 1.357674400 4.073023100 4.073023100
      # Sb 4.073023100 4.073023100 1.357674400
    # &END COORD
    &TOPOLOGY
      &CENTER_COORDINATES
  	  &END CENTER_COORDINATES
 	  COORD_FILE_NAME sisb.cif
  	  COORD_FILE_FORMAT CIF # 除了 cif，还可以读入 pdb、xtl，xyz 等格式；
    &END TOPOLOGY
    &KIND Sb
      ELEMENT Sb
      BASIS_SET TZVP-MOLOPT-SR-GTH-q5
      POTENTIAL GTH-PBE
    &END KIND
    &KIND Si
      ELEMENT Si
      BASIS_SET TZVP-MOLOPT-GTH-q4
      POTENTIAL GTH-PBE
    &END KIND
  &END SUBSYS
  &DFT
	@SET DATAPATH /home/changruiwang-ICME/Software/cp2k-2023.1/data
    BASIS_SET_FILE_NAME ${DATAPATH}/BASIS_MOLOPT # 基组文件 1；
    BASIS_SET_FILE_NAME ${DATAPATH}/BASIS_MOLOPT_UCL # 基组文件 2，相比 BASIS_MOLOPT，BASIS_MOLOPT_UCL 可选的元素种类更多；
    POTENTIAL_FILE_NAME ${DATAPATH}/POTENTIAL # 赝势文件；
    # WFN_RESTART_FILE_NAME cp2k-RESTART.wfn # 读取波函数；
    # CHARGE 0 # 体系整体电荷；
    # MULTIPLICITY 1 # 体系的整体自旋多重度；
    &QS
      EPS_DEFAULT 1.0E-14
      EXTRAPOLATION ASPC
      EXTRAPOLATION_ORDER 3
    &END QS
    &POISSON
      PERIODIC XYZ
 	  PSOLVER PERIODIC
    &END POISSON
    &XC
      &XC_FUNCTIONAL PBE # 控制泛函；
      &END XC_FUNCTIONAL
      &VDW_POTENTIAL # 控制色散校正；
        POTENTIAL_TYPE PAIR_POTENTIAL # 除了 PAIR_POTENTIAL 还有 NON_LOCAL（例如 RVV10）；
   	    &PAIR_POTENTIAL
   		  TYPE DFTD3(BJ)
          PARAMETER_FILE_NAME ${DATAPATH}/dftd3.dat
      	  REFERENCE_FUNCTIONAL PBE
		  R_CUTOFF 15 # 截断半径；
      	  # CALCULATE_C9_TERM T # C9 项，同时增加计算精度和计算量；
		  &PRINT_DFTD
	        FILENAME cp2k-dftd3.out
		  &END PRINT_DFTD
 	    &END PAIR_POTENTIAL
      &END VDW_POTENTIAL
    &END XC
    &MGRID
      NGRIDS 5 # 对 MOLOPT-GTH 基组，5 是最优设置；
      CUTOFF 400 # 取决于元素种类；
      REL_CUTOFF 60 # 默认 40，一般取 50 或 60；
      # USE_FINER_GRID T # 用于提高网格精细度；
    &END MGRID
    &SCF
      MAX_SCF 300 # SCF 迭代上限，类似于 VASP 中的 NELM；
      EPS_SCF 1.0E-6 # SCF 能量收敛标准，默认是 1E-5，单位是 hartree；
      SCF_GUESS RESTART
      &DIAGONALIZATION
        ALGORITHM STANDARD
        EPS_ADAPT 0.01
      &END DIAGONALIZATION
      &MIXING
        METHOD BROYDEN_MIXING
        ALPHA 0.4
        BETA 1.5
        NBROYDEN 8
      &END MIXING
      &SMEAR
        METHOD FERMI_DIRAC
  	    ELECTRONIC_TEMPERATURE 300
      &END SMEAR
      ADDED_MOS 500 # 求解一些额外的空轨道从而能被热激发的电子所占据；
      &PRINT
        &RESTART
          FILENAME cp2k-RESTART.wfn
   		  BACKUP_COPIES 0 
  	    &END RESTART
      &END PRINT
    &END SCF
    &PRINT
      &E_DENSITY_CUBE
        FILENAME cube
     	STRIDE 1 1 1 #Stride of exported cube file
   	  &END E_DENSITY_CUBE
    &END PRINT
  &END DFT
  &PRINT 
    &FORCES ON
      FILENAME cp2k-FORCE.dat
    &END FORCES 
  &END PRINT
&END FORCE_EVAL
```
CP2K 通常默认使用单 Gamma 点计算，需要使用足够大的晶胞以减少周期性误差。
如果晶胞太小，部分基组函数可能超出晶胞边界，导致重叠矩阵求逆过程出现数值问题，从而使得结果不可靠。
可以考虑增大晶胞尺寸或在输入文件中调整边界条件（如 `POISSON` 模块中的 `PERIODIC NONE`），
并适当提高平面波基组的截断能量（`CUTOFF`）和相对截断能量（`REL_CUTOFF`）。

切记，不能直接使用 VASP 中的小晶胞来进行 CP2K 的计算。
需要使用多 k 点时，可以在 `KPOINTS` 中指定网格的 k 点采样。

CP2K 计算中需要指定赝势和基组文件（类似于 VASP 中的 `POTCAR`）。
赝势文件通常选择位于 `GTH_POTENTIALS` 目录下的 `PBE` 或 `BLYP` 泛函文件，并根据元素的价电子数选择不同的 `-q` 后缀，如 `-q4` 表示该赝势包含 4 个价电子。
基组文件常用 `BASIS_MOLOPT`，它经过专门的优化，适用于大部分元素体系；
而 `BASIS_MOLOPT_UCL` 是伦敦大学学院提供的扩展基组，适用于更多类型的元素。

在基组文件中，`SZV`（Single Zeta Valence）、`DZVP`（Double Zeta Valence with Polarization）、`TZVP`（Triple Zeta Valence with Polarization）分别代表基组的劈裂（多少层电子壳层）和极化函数的添加情况。
其中，`VP`（Valence Polarized）表示基组中包含极化函数。通常，劈裂层数和极化函数越多，计算结果越精确，但计算量也会显著增加。
计算量的顺序一般为：`SZV < DZVP < TZVP < TZV2P < TZV2PX`。

选择基组时应注意与赝势中的价电子数保持一致。
基组的大小对基组重叠误差（Basis Set Superposition Error, BSSE）有显著影响。
较大的基组（如 `TZVP` 或 `TZV2P`）通常可以有效降低 BSSE 误差，因此对于大体系或精度要求较高的计算，推荐使用 `DZVP` 或 `TZVP`。
在计算过程中，可以使用 COUNTERPOISE 方法来校正 BSSE 误差，确保结果的精度。

---

## FORCE_EVAL 中的 SCF 收敛算法

CP2K 提供了两种常用的 SCF（Self-Consistent Field）收敛算法：
1. 对角化算法（`DIAG`）：基于直接对角化的方式求解 KS（Kohn-Sham）方程。对于体系的带隙很小（如半导体或金属）或几乎没有（如金属体系），推荐使用 `DIAG` 算法，并开启电子展宽（`SMEAR`）来帮助收敛。

注意：在使用 `DIAG` 算法时，必须设置 `ADDED_MOS` 参数来指定额外的虚拟轨道数量，以防止出现电子填充不完全的问题，并提升收敛效果。
`ADDED_MOS` 的值通常需要设置为体系中价带电子数的 `5-10%` 左右，但对于金属体系或需要处理激发态的计算，可能需要更多的虚拟轨道以确保所有电子可以正确填充。

2. 轨道变换算法（`OT`）：一种基于轨道变换的 SCF 优化方法，通常比对角化算法更高效，尤其适用于具有较大带隙的绝缘体或介电材料体系。在体系中若带隙较大或者体系较为复杂时，`OT` 算法往往能显著加快 SCF 的收敛。

在使用 `DIAG` 或 `OT` 算法时，设置合适的混合参数（`MIXING`）是加速收敛的关键。
通常，对于金属体系（小带隙或无带隙），较低的混合参数（如 `MIXING 0.1`）可以避免 SCF 振荡；
而对于带隙较大的绝缘体，较高的混合参数（如 `MIXING 0.4`）则有助于加快收敛。

`DIAG` 算法的输入文件参考上一节。
`OT` 算法的 `.inp` 输入文件如下所示（只展示 `FORCE_EVAL/DFT` 下的 `SCF`），替换上一节中 `SCF` 即可：
```
    &SCF
      MAX_SCF 30
      EPS_SCF 1E-06
      SCF_GUESS RESTART
      &OT
	    ALGORITHM IRAC
        MINIMIZER CG
        LINESEARCH 3PNT
        PRECONDITIONER FULL_SINGLE_INVERSE
      &END OT
      &OUTER_SCF
        EPS_SCF 1.0E-05
        MAX_SCF 5
      &END OUTER_SCF
      &PRINT
        &RESTART OFF
          BACKUP_COPIES 0 
        &END RESTART
      &END PRINT
    &END SCF
```

`OT/MINIMIZER` 常用的设置有 `CG`、`DIIS` 以及 `BROYDEN`。
其中，`CG` 是最为稳定的算法，一般的计算任务都可以使用 `CG` 算法，并设置 `OT/LINESEARCH` 为比默认 `2PNT` 更贵但也更稳健的 `3PNT`。
`DIIS` 算法速度比较快，但没有 `CG` 稳定。
如果 `CG` 算法和 `DIIS` 算法收敛都有问题时，可以尝试使用 `BROYDEN` 算法。
采用 `OT` 时，推荐开启 `SCF/OUTER_SCF`，这是加速收敛的一种方法。
开启 `OUTER_SCF` 时 `SCF/MAX_SCF` 应当非常小，`15` 至 `35` 是比较合适的范围。
`OUTER_SCF` 迭代圈数通过 `SCF/OUTER_SCF/MAX_SCF` 控制，一般设为 `5` 即可。
实际 SCF 迭代次数上限等于 `SCF/MAX_SCF *（1+SCF/OUTER_SCF/MAX_SCF）`。
另外，`SCF/EPS_SCF` 设置的是 SCF 总的收敛标准（一般是 `1E-6` ），而 `SCF/OUTER_SCF/EPS_SCF` 设置的是 `OUTER_SCF` 的收敛标准（一般是 `1E-5`），后者绝对值应当大于或等于前者。
`OT/PRECONDITIONER` 中，`FULL_ALL` 通常最稳定，但耗时最长（`GAPW` 计算只能使用 `FULL_ALL` 另说）；
大体系可以尝试 `FULL_SINGLE_INVERSE` 和 `FULL_KINETIC`。
`OT/ALGORITHM` 可以从默认的 `STRICT` 改为 `IRAC`，SCF 收敛会更稳健。
但实际上，即使是对于非金属体系，有时候 `DIAG` 算法也可能比 `OT` 算法速度更快。
所以，在进行大规模的计算之前最好进行充分的测试。

---
