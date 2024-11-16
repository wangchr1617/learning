
# QE 结构优化

下面以 GeTe 为例，讲述如何使用 QE 优化晶格参数和原子坐标。

本例从 GeTe.cif 文件出发。使用 `Multiwfn_noGUI GeTe.cif` 打开 Multiwfn 程序，依次键入 `QE`、`pw.inp`、`q` 生成一个基础的输入文件模板（内含结构信息）如下所示：

```
&control
 calculation= 'scf'
 prefix= 'GeTe'
 pseudo_dir= './'
/
&system
 ibrav= 0
 nat=     6
 ntyp=    2
 ecutwfc= 44.0
/
&electrons
/
CELL_PARAMETERS angstrom
    4.21855598    0.00000000    0.00000000
   -2.10927799    3.65337664    0.00000000
    0.00000000    0.00000000   11.03915739
ATOMIC_SPECIES
Ge  72.6276  Pseudopotential_file
Te 127.6031  Pseudopotential_file
ATOMIC_POSITIONS angstrom
Ge    0.000000    0.000000   10.968286
Ge    2.109299    1.217780    3.608811
Ge   -0.000021    2.435597    7.288493
Te    2.109299    1.217780    9.502065
Te   -0.000021    2.435597    2.142590
Te    0.000000    0.000000    5.822272
K_POINTS automatic
3 3 3 0 0 0
``` 

然后在这个模板的基础上增加结构优化相关的 section 如下所示：

```
&CONTROL
 calculation = 'vc-relax' 	! A string describing the task to be performed.
 verbosity = 'high'
 nstep = 100 				! number of molecular-dynamics or structural optimization steps performed in this run.
 tstress=.true. 			! calculate stress.
 tprnfor=.true. 			! calculate forces.
 outdir= './tmp' 			! input, temporary, output files are found in this directory
 prefix = 'pwscf' 			! prepended to input/output filenames
 forc_conv_thr = 1.0d-4 	! Convergence threshold on forces (a.u) for ionic minimization
 disk_io = 'low' 			! Specifies the amount of disk I/O activity
 pseudo_dir= './' 			! directory containing pseudopotential files
/
&SYSTEM
 ibrav = 0 					! Bravais-lattice index.
 nat = 6 					! number of atoms in the unit cell
 ntyp = 2 					! number of types of atoms in the unit cell
 ecutwfc= 50 				! kinetic energy cutoff (Ry) for wavefunctions
 ecutrho = 500 				! Kinetic energy cutoff (Ry) for charge density and potential
 nosym = .false. 			! if (.TRUE.) symmetry is not used.
 noinv = .false. 			! if (.TRUE.) disable the usage of k => -k symmetry
 occupations = 'smearing'
 degauss = 1.0d-9 			! value of the gaussian spreading (Ry) for brillouin-zone
 smearing = 'gauss'
 nspin = 1 					! non-polarized calculation
 vdw_corr = 'dft-d3' 		! Type of Van der Waals correction
 dftd3_version = 4 			! Grimme-D3 BJ damping
 /
&ELECTRONS
 electron_maxstep = 100 	! maximum number of iterations in a scf step.
 scf_must_converge = .false. ! If .false. do not stop molecular dynamics or ionic relaxation when electron_maxstep is reached.
 conv_thr = 1.0d-7 			! Convergence threshold for selfconsistency
 mixing_mode = 'plain' 		! charge density Broyden mixing
 mixing_beta = 0.8d0 		! mixing factor for self-consistency
 diagonalization = 'david' 	! Davidson iterative diagonalization with overlap matrix
 tqr = .false. 				! If .true., use a real-space algorithm for augmentation charges of ultrasoft pseudopotentials and PAWsets.
 real_space = .false. 		! If .true., exploit real-space localization to compute matrix elements for nonlocal projectors.
/
&IONS
 ion_dynamics = 'bfgs' 		! use BFGS quasi-newton algorithm
 trust_radius_min = 1.d-3   ! Minimum ionic displacement
 trust_radius_ini = 0.5d0   ! Initial ionic displacement in the structural relaxation
/
&CELL
 cell_dynamics = 'bfgs' 	! use BFGS quasi-newton algorithm
 press = 0.0 				! Target pressure [KBar] in a variable-cell md or relaxation run.
 press_conv_thr = 0.1 		! Convergence threshold on the pressure for variable cell
 cell_dofree = 'all' 		! all axis and angles are moved
/
CELL_PARAMETERS angstrom
    4.21855598    0.00000000    0.00000000
   -2.10927799    3.65337664    0.00000000
    0.00000000    0.00000000   11.03915739
ATOMIC_SPECIES
Ge  72.6276  Ge.UPF
Te 127.6031  Te.UPF
ATOMIC_POSITIONS angstrom
Ge    0.000000    0.000000   10.968286
Ge    2.109299    1.217780    3.608811
Ge   -0.000021    2.435597    7.288493
Te    2.109299    1.217780    9.502065
Te   -0.000021    2.435597    2.142590
Te    0.000000    0.000000    5.822272
K_POINTS automatic
12 12 3 0 0 0
```

其中赝势文件 `Ge.UPF` 和 `Te.UPF` 是 SSSP 赝势库中拷贝到本地的，`K_POINTS` 对应 VASPKIT 设置 Accuracy Levels 为 0.03 时的 KPOINTS 文件。

使用命令 `qsub runqe.pbs` 提交作业。脚本如下所示：

```
#PBS -N qe
#PBS -l nodes=1:ppn=32
#PBS -l walltime=24:00:00
#PBS -q manycores
#PBS -V
#PBS -S /bin/bash

export MODULEPATH=/opt/modulefiles
module load intel/2020.1.217
module load gcc/9.3

NP=`cat $PBS_NODEFILE | wc -l`
NN=`cat $PBS_NODEFILE | sort | uniq | tee /tmp/nodes.$$ | wc -l`

cat $PBS_NODEFILE > /tmp/nodefile.$$
cd $PBS_O_WORKDIR
ulimit -s unlimited

EXEC=/home/changruiwang-ICME/Software/qe-7.3.1/bin/pw.x
mpirun -machinefile $PBS_NODEFILE -np $NP $EXEC -pd .true. < pw.inp > output
```

---

## 结果分析与讨论

计算过程中使用下列命令可以查看优化过程中的收敛情况：

```
grep " Total force =" output # 查看优化过程中各原子总的受力情况
grep -A 20 "Forces acting on atoms (cartesian axes, Ry/au):" output # 查看优化过程中各原子的受力情况
grep ! output # 查看优化过程中能量的变化过程
grep -A 12 " Total force =" output # 查看优化过程中压力张量大小
```

计算结束后，运行

```
awk '/Begin final coordinates/,/End final coordinates/{print $0}' output
```

得到 vc-relax 的结果如下所示：

```
Begin final coordinates
     new unit-cell volume =   1148.92721 a.u.^3 (   170.25345 Ang^3 )
     density =      5.85875 g/cm^3

CELL_PARAMETERS (angstrom)
   4.219080033  -0.000000049  -0.000000162
  -2.109540059   3.653830459   0.000000162
  -0.000000424   0.000000245  11.044085481

ATOMIC_POSITIONS (angstrom)
Ge               0.0000043686       -0.0000025222       10.9724916280
Ge               2.1095476175        1.2179389741        3.6100565748
Ge              -0.0000125968        2.4358946620        7.2916713080
Te               2.1095460691        1.2179398681        9.5065439808
Te              -0.0000001133        2.4358874547        2.1438972863
Te              -0.0000068854        0.0000039753        5.8254150101
End final coordinates
```
