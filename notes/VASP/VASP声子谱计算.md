
# VASP 声子谱计算

## Phonopy 安装

`Phonopy` 是计算声子谱的主流工具。推荐通过以下方法安装：

### 在线安装
使用 `conda install -c conda-forge phonopy` 安装。

### 离线安装
考虑到集群无法联网的情况，可以通过以下步骤离线安装：

1. **下载 Phonopy 源码或安装包**: 
   从 [PyPI](https://pypi.org/project/phonopy/) 或 [GitHub](https://github.com/phonopy/phonopy.git) 获取最新版本的安装包，并解压。
   
2. **安装**:
   ```
   cd phonopy
   pip install .
   ```

3. **安装依赖**:
   如果遇到依赖问题，可逐一手动安装。
   以 `Phonopy` 的依赖库 `spglib` 为例，可以访问 [PyPI](https://pypi.org/project/spglib/) 下载目标版本的 `.whl` 文件，并安装：
   ```
   pip install spglib-2.0.2-cp37-cp37m-manylinux1_x86_64.whl
   ```
   注意根据集群修改文件名以符合当前集群操作系统规范。
   例如，`spglib-2.0.2-cp37-cp37m-manylinux_2_17_x86_64.manylinux2014_x86_64.whl` 需要改为 `spglib-2.0.2-cp37-cp37m-linux_x86_64.whl`才能在集群 2 安装。

### 集群兼容性注意事项
- **低版本 Anaconda**: 对于旧版本集群（如集群 2），推荐安装兼容的 `Phonopy` 版本（如 2.14.0）。
- **GLIBC 限制**: 若集群支持的 GLIBC 版本较低（如仅支持到 2.12，可通过 `strings /lib64/libc.so.6 | grep GLIBC` 查看），确保使用兼容版本的依赖库（如 `spglib-1.16.0`）。

### 验证安装
安装完成后，运行 `phonopy`，若成功，输出类似以下信息：
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

若显示 `bash: command not found`，需将 `Phonopy` 的路径添加到环境变量：
```
echo 'export PATH=/your_phonopy_path/bin:$PATH' >> ~/.bashrc
source ~/.bashrc
```

---

## 声子谱计算

VASP 结合 `Phonopy` 计算声子谱有两种主要方法：

1. **有限位移法**: 又称 frozen-phonon 方法或直接法。

> 有限位移法是通过在优化后的平衡结构中引入原子位移，计算作用在原子上的 Hellmann-Feynman 力，进而由动力学矩阵算出声子色散曲线。
> 有限位移法的缺陷在于它要求声子波矢与原胞边界 supersize 正交，或者原胞足够大使得 Hellmann-Feynman 力在原胞外可以忽略不计。
> 这使得对于复杂系统，如对称性高的晶体、合金、超晶格等材料需要采用超原胞。
> 超原胞的采用使计算量急剧增加，极大的限制了该方法的使用。

2. **密度泛函微扰理论 (DFPT)**。

> 1987年，Baroni、Giannozzi 和 Testa 提出了 DFPT 方法。
> DFPT 通过计算系统能量对外场微扰的响应来求出晶格动力学性质。
> 该方法最大的优势在于它不限定微扰的波矢与原胞边界正交，不需要超原胞也可以对任意波矢求解。
> Castep、Vasp 等采用的是一种 linear response theory 的方法（或者称为 density perturbation functional theory，DFPT）。
> 直接计算出原子的移动而导致的势场变化，再进一步构造出动力学矩阵，进而计算出声子谱。

无论使用哪种方法，都需要非常高的计算精度（如 `EDIFF = 1E-08`），并确保初始结构已充分优化以减少虚频的可能性。

### Phonopy 前处理

- 声子谱计算要求超胞的晶格常数通常大于 10 Å，原子数在 100 左右。
- 优化完成后，执行扩胞操作，将 POSCAR 沿着 a、b、c 轴分别扩 4 倍：
  ```
  phonopy -d --dim="4 4 4" --pa="AUTO"
  ```
  生成以下文件：
  - `SPOSCAR`: 完美超胞的结构文件。
  - `phonopy_disp.yaml`: 原子位移信息。
  - `POSCAR-{number}`: 包含位移的结构文件。

### 有限位移法

1. 使用脚本将位移结构文件分配到不同的文件夹进行单点能计算：
   ```bash
   # disp.sh
   for i in POSCAR-0*
   do
     name=$(echo "$i" | grep -oE '[0-9]+')
     mkdir $name
     mv $i $name/POSCAR
   done
   ```
2. 在各个文件夹中运行高精度单点计算（`IBRION = -1`），生成 `vasprun.xml`。
3. 使用命令 `phonopy -f {001..004}/vasprun.xml --pa="AUTO"` 生成 `FORCE_SETS` 文件，`--pa` 参数适用于惯用胞计算；
   或者使用命令 `phonopy --fc {001..004}/vasprun.xml` 生成 `FORCE_CONSTANTS` 文件。

### 密度泛函微扰理论 (DFPT)

1. 将 `SPOSCAR` 复制为 `POSCAR`。
2. 在 `INCAR` 中设置 `IBRION = 8`、`NSW = 1`、`ISIF = 3` 并注释掉并行参数。若需要考虑色散校正或自旋轨道耦合，可设置 `IBRION = 5` 或 `IBRION = 6`。
3. 使用命令 `phonopy --fc vasprun.xml` 生成 `FORCE_CONSTANTS` 文件。

### 检查输出

1. **检查计算**:
   ```
   grep hessian vasprun.xml
   grep THz OUTCAR # 查看虚频
   ```
2. **生成 FORCE_CONSTANTS**:
   如果需要从 FORCE_SETS 获得 FORCE_CONSTANTS。创建 `writefc.conf` 文件：
   ```
   FORCE_CONSTANTS = WRITE
   DIM = 4 4 4
   ```
   运行命令：
   ```
   phonopy writefc.conf
   ```

### Phonopy 后处理绘制声子谱

1. **生成 `band.conf`**:
   使用 `vaspkit 305 3` 创建 `KPATH.phonopy` 文件，并修改为 `band.conf`。
   然后对 `band.conf` 文件稍作修改，例如：
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
   FC_SYMMETRY = .TRUE.
   
   # FORCE_SETS = READ # 也可以选择读取 FORCE_SETS
   # IRREPS = 0  0  0
   # SHOW_IRREPS = .TRUE.
   # LITTLE_COGROUP = .TRUE.
   ```
   
   `band.conf` 文件的参数含义详见 Phonopy 官网。
   
2. **运行绘图**:
   依次运行下面命令输出 `band.yaml`、`band.dat` 和 `band_dos.pdf`：
   ```
   phonopy --dim="4 4 4" -c POSCAR --factor=521.471 -p -s band.conf
   phonopy-bandplot --gnuplot band.yaml > band.dat
   ```
   其中，`--factor=521.471` 意指将单位从 THZ 换为 $cm^{-1}$ 。
   
   `band.dat` 可以直接拖进 Origin 画图，第 2 行是 band.conf 中高对称点对应的横坐标。
   访问 [网址](https://henriquemiranda.github.io/phononwebsite/phonon.html) 并上传 `band.yaml` 可以查看声子谱的振动模式。
   