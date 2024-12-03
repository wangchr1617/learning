
# LAMMPS 计算声子谱

## 安装环境

首先 `pip install phonolammps` 安装 phonolammps 。
然后安装 LAMMPS 的 Python 接口，依次运行：
```
make serial mode=shlib
make install-python
```

在 Python 能 `from lammps import lammps` 导入即可。

---

## 使用 NEP 势函数计算 GeTe 声子谱

然后准备 LAMMPS 的输入文件 in.lammps 如下所示：

```
units    metal
atom_style  atomic
dimension 3
boundary p p p
read_data     graphene.data #读取结构
pair_style    nep nep.txt  #选择自己的势函数
pair_coeff    * *
```

利用 phonolammps 生成二阶力矩阵：

```
phonolammps in.lammps -c POSCAR --dim 4 4 4
```

运行完之后会生成二阶力常数矩阵 `FORCE_CONSTANTS`。

然后利用 `vaspkit 305 3` 生成高对称点路径文件 `KPATH.in` 并重命名为 `band.conf`。
修改 `band.conf` 如下所示：
```
ATOM_NAME = Ge Te # 元素名称
NPOINTS = 501
DIM = 4 4 4 # 扩胞指数
BAND = 0.000000 0.000000 0.000000 0.500000 0.500000 0.500000 0.759157 0.240843 0.500000, 0.500000 -0.240843 0.240843 0.500000 0.000000 0.000000 0.000000 0.000000 0.000000 0.370421 -0.370421 0.000000, 0.629579 0.000000 0.370421 0.500000 0.000000 0.500000 0.000000 0.000000 0.000000
BAND_LABELS = $\Gamma$ T H$_2$ H$_0$ L $\Gamma$ S$_0$ S$_2$ F $\Gamma$

MP = 30 30 30
TETRAHEDRON = .TRUE.
PDOS = 1, 2
BAND_CONNECTION = .TRUE.
EIGENVECTORS = .TRUE. # 输出本征矢
FORCE_CONSTANTS = READ # 读取 FORCE_CONSTANTS

# FORCE_SETS = READ # 也可以选择读取 FORCE_SETS
# IRREPS = 0  0  0
# SHOW_IRREPS = .TRUE.
# LITTLE_COGROUP = .TRUE.
```

依次运行下面命令输出 `band.yaml`、`band.dat` 和 `band_dos.pdf`：
```
phonopy --dim="4 4 4" -c POSCAR --factor=521.471 -p -s band.conf
phonopy-bandplot --gnuplot band.yaml > band.dat
```
其中，`--factor=521.471` 意指将单位从 THZ 换为 $cm^{-1}$ 。

`band.dat` 可以直接拖进 Origin 画图，第 2 行是 band.conf 中高对称点对应的横坐标。
访问 [网址](https://henriquemiranda.github.io/phononwebsite/phonon.html) 并上传 `band.yaml` 可以查看声子谱的振动模式。
