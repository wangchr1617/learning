# CP2K_入门

CP2K 的输入文件一般命名为 `cp2k.inp`，其中关键词以 section 和 subsection 的形式，一环套一环，如下图所示。

<div align="left">
<img src="./figures/fig_001.png" width = "25%" />
</div>

每一个 section 都是以 `&section` 开头，以 `&END section` 结尾，顺序随意，但是嵌套不能乱。
在 section 中，每行只写一个关键词，关键词后面接参数，CP2K对大小写和空格不敏感。

常用的 section 包括 `GLOBAL`、`FORCE_EVAL`、`MOTION` 等，下面将逐个介绍它们。

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
5. `CELL_OPT`，晶胞优化，常与GEO_OPT联用；
6. `BAND`，NEB计算；

`PRINT_LEVEL` 控制输出详略，一般设置成 `LOW` 或 `MEDIUM` 即可。

CP2K 中计算能量需要设置 `GLOBAL/RUN_TYPE` 为 `ENERGY`；
如果还要计算受力，需要设置 `GLOBAL/RUN_TYPE` 为 `ENERGY_FORCE`，
并设置 `FORCE_EVAL/PRINT/FORCES` 为 `ON`。
如果要将受力信息存储到文件中，则需要设定 `FROCE_EVAL/PRINT/FORCES/FILENAME` 文件名；
否则受力信息将会打印到 `*.out` 文件中。

注意，在 CP2K 中，默认的能量和距离单位分别是 Hartree (Hartree 原子单位) 和 Bohr (波尔半径)。
如果需要将它们转换成常用的 eV（电子伏特）和 Å（埃）单位，可以使用以下换算关系：
```
1 Hartree = 27.2114 eV
1 Bohr = 0.529177 Å
```
反之：
```
1 eV = 0.0367493 Hartree
1 Å = 1.88973 Bohr
```
