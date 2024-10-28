# VASP 批量计算

本节介绍几个常用于 VASP 批量计算的脚本：`POSCAR.sh`、`link.sh`、`batchrunvaspkit.sh`、`batchrunvasp.pbs`、`check.sh`、`find.sh` 和 `CONTCAR.sh`，
并结合示例讲述如何构建 `.xyz` 格式的势函数训练集。

## 生成结构

准备训练集之前，首先要明确做一个怎样的势函数？
通用势函数需要采样非常丰富的结构；
而专注于研究某种性质的势函数，采样覆盖目标问题相关的势能面区域即可。

常规的势函数训练集包括对应温度区间的 AIMD 采样结构，结合一定的原子坐标微扰结构即可。
AIMD 的精度不用要求太高，所以如果 AIMD 的温度区间是低温、室温或中温，使用 VASP 构建训练集时可以使用 MLFF（VASP 版本大于等于 6.3）加速采样。

模拟过程中，AIMD 一般间隔 `0.1 ps` 对轨迹抽样，更小时间间隔内的结构之间具有强关联性，对增加训练集多样性无益，徒增训练集大小。

如果研究高压或辐照体系，采样之后的单点能计算可能需要非常高的精度（例如 `KSPACING = 0.1`、`SIGMA = 0.02`）来确保原子受力的准确。
此外，辐照相关的势函数还需要非常多的小团簇结构，以充分覆盖原子局域环境的剧烈变化。

### VASP MLFF 加速采样

INCAR 设置如下所示：
```
NCORE = 12
ISTART = 0
ICHARG = 2
LWAVE = .F.
LCHARG = .F.
KSPACING = 2 # 使用较大的 KSPACING 加速模拟
KGAMMA = .T.
ISYM = 0
IVDW = 12 # 根据需要启用色散校正

ALGO = Normal # Fast 更快，但是 Normal 更稳定
NELM = 120
NELMIN = 4
EDIFF = 1E-05 # 使用较低的收敛标准加速模拟

IBRION = 0
ISIF = 2
POTIM = 2.0 # 设置模拟步长
NSW = 5000 # 设置模拟步数
TEBEG = 300 # 设置模拟温度
TEEND = 300 # 设置模拟温度
MDALGO = 2
SMASS = 0
NBLOCK = 1

# MLFF 相关，高温采样时全部注释即可
ML_MODE = train
ML_LMLFF = .T.
# ML_ISTART = 1
# ML_CX = -0.1 # 根据需要增大第一性原理计算比率
# ML_IALGO_LINREG = 3
# ML_RCUT2 = 6.0
# ML_RCUT1 = 6.0
# ML_CTIFOR = 1000
# ML_MB = 3000 # 模拟失败时根据LCONF适当增大
# ML_LBASIS_DISCARD = .TRUE.

ENCUT = 400 # 减小 ENCUT 加速模拟
ISMEAR = 0
SIGMA = 0.05
PREC = N
LREAL = A
ADDGRID = .T.
```

设置 MLFF 时，使用命令 `grep STATUS ML_LOGFILE` 查看计算进程，命令 `tail -f ML_LOGFILE | grep STATUS` 实时监测进程。
计算完成后，使用 `grep ERR ML_LOGFILE > mlff_err` 输出整个过程的误差变化；
使用 `grep ‘free  energy ML TOTEN’ OUTCAR > mlff_ene` 得到整个过程的能量变化。
如果想要续算 MLFF，在当前目录下
```
cp ML_ABN ML_AB
cp ML_FFN ML_FF
cp CONTCAR POSCAR
```
并在 INCAR 中设置 `ML_ISTART = 1` 即可。

---