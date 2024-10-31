# VASP 声子谱计算

贡献者：胡京京、王昌锐

---

## Phonopy 安装

`Phonopy` 是计算声子谱主流的软件。官网推荐的安装方式是使用 `conda install -c conda-forge phonopy` 安装，
但是考虑到课题组集群不能联网的现实条件，本节首先介绍如何在自己的账号下离线安装 `Phonopy`。

1. 下载 Phonopy 包

访问 [PyPI](https://pypi.org/project/phonopy/) 或者 [Github](https://github.com/phonopy/phonopy.git) 下载 `Phonopy` 安装包并解压。

2. 安装 Phonopy

进入下载的目录：
```
cd phonopy
```
 
然后运行：
```
pip install .
```
如果自动解析并安装依赖失败，可以考虑手动安装所需依赖并通过 `python setup.py install` 直接调用 `setuptools`，使用 `setup.py` 中定义的配置来安装。

3. 离线安装依赖库

离线安装软件的原则是 “缺啥补啥” 。
以 `Phonopy` 的依赖库 `spglib` 为例，可以访问 [PyPI](https://pypi.org/project/spglib/) 下载目标版本的 `.whl` 文件。

注意根据集群修改文件名以符合当前集群操作系统规范。
例如，`spglib-2.0.2-cp37-cp37m-manylinux_2_17_x86_64.manylinux2014_x86_64.whl` 需要改为 `spglib-2.0.2-cp37-cp37m-linux_x86_64.whl`才能使用以下命令安装：
```
pip install spglib-2.0.2-cp37-cp37m-linux_x86_64.whl
```

4. 针对集群 2 的兼容性注意事项

集群 2 因为支持的 Anaconda 版本过低，不能下载最新版本的 `Phonopy`，确定 `2.14.0` 版是兼容的。
另外，由于集群 2 的 `GLIBC` 版本限制（仅支持到 2.12，可通过 `strings /lib64/libc.so.6 | grep GLIBC` 查看当前支持的 `GLIBC` 版本），
`spglib` 也不能安装最新版本，目前确定 `spglib-1.16.0-cp37-cp37m-manylinux1_x86_64.whl` 是兼容的。

5. 测试安装

安装完成后，在命令行输入 `phonopy`。如果安装成功，会显示如下界面，标注 `Python` 和 `Spglib` 的版本号：
```
        _
  _ __ | |__   ___  _ __   ___   _ __  _   _
 | '_ \| '_ \ / _ \| '_ \ / _ \ | '_ \| | | |
 | |_) | | | | (_) | | | | (_) || |_) | |_| |
 | .__/|_| |_|\___/|_| |_|\___(_) .__/ \__, |
 |_|                            |_|    |___/
                                      2.12.0

Python version 3.9.6
Spglib version 1.16.2
```

6. 配置环境变量

如果输入 `phonopy` 时出现 `bash: command not found` 的提示，
则需要将 `Phonopy` 可执行文件所在的 `bin` 文件夹路径手动添加到 `~/.bashrc` 的环境变量里，
方式如下：
```
echo 'export PATH=/your_phonopy_path/bin:$PATH' >> ~/.bashrc
source ~/.bashrc
```

---

## 声子谱计算

1. 准备进行过结构优化的 POSCAR
2. `phonopy -d --dim="2 2 1" ` # 在 xyz 方向扩胞大小
3. 此时应该有：phonopy_disp.yaml  POSCAR-001  POSCAR-003  POSCAR-005  POSCAR-002  POSCAR-004  POSCAR-006  SPOSCAR POSCAR 这些文件，然后

```cp POSCAR POSCAR_unit 
  mv SPOSCAR POSCAR # 生成了 SPOSCAR文件，把原来的 POSCAR 保存为 POSCAR_unit，把SPOSCAR 改名为 POSCAR 用于后续计算
```
4. 准备 INCAR

```
IOPTCELL = 1 1 0 1 1 0 0 0 0 # 控制z方向不变
ADDGRID = True
ALGO = Normal
EDIFF = 1e-08
ENCUT = 600
IBRION = 8
ISIF = 3
ISMEAR = 0
ISPIN = 1
ISYM = 2
LASPH = True
LCHARG = False
LORBIT = 11
LREAL = False
LWAVE = False
NELM = 100
NSW = 1
POTIM = 0.01
PREC = Normal
SIGMA = 0.1
EDIFFG = -1e-03
```
5. 生成KPOINTS，方法和做结构优化的一样，复制pbs文件，提交
6. 算完之后，vaspkit 305 2 —— 生成KPATH.phonopy 和 band.conf
7. 
```grep hession vasprun.xml # 显示 说明计算成功
phonopy --fc vasprun.xml # 获得力常熟文件
grep THz OUTCAR # 查看虚频
cp KPATH.phonopy band.conf 
```

8. `vi band.conf` 该文件如下
```
ATOM_NAME = B Cr Mo # 这里原来没有，需要加上原子种类
NPOINTS = 501
DIM = 2 2 1 #在 xyz 方向扩胞大小
BAND = 0.000000 0.000000 0.000000 0.500000 0.000000 0.000000 0.333333 0.333333 0.000000 0.000000 0.000000 0.000000
BAND_LABELS = $\Gamma$ M K $\Gamma$

MP = 30 30 30
TETRAHEDRON = .TRUE.
#PDOS = 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50, 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75, 76 77 78 79 80 81 82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99 100 #是否计算PDOS，算的话取消注释
BAND_CONNECTION = .TRUE.
FORCE_CONSTANTS = READ

# FORCE_SETS = READ
# IRREPS = 0  0  0
# SHOW_IRREPS = .TRUE.
# LITTLE_COGROUP = .TRUE.
```
9. 
```phonopy --dim="2 2 1" -c POSCAR_unit band.conf
mv POSCAR SPOSCAR
mv POSCAR_unit POSCAR
phonopy --factor=521.471 band.conf
phonopy -bandplot --gnuplot band.yaml > band.dat
cp ~/script/phonbydat.py .
python phonbydat.py band.dat
```
10. 绘制声子谱如图
<div align="left">
<img src="./figures/002.png" width = "50%" />
</div>

