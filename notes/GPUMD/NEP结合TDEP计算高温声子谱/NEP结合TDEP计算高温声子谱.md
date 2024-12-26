
# TDEP 计算高温声子谱

本教程将介绍 TDEP（Temperature Dependent Effective Potential）方法的基础知识。
内容包括如何从分子动力学模拟生成的位置-力数据中提取有效的原子间力常数（IFCs），以及如何进一步计算声子色散关系、态密度和振动自由能。

所用数据来自不同温度下、使用包含 128 个原子的菱方相 GeTe 4 × 4 × 4 超胞进行的分子动力学模拟。
首先在 NPT 系综下使用超大胞（ > 10000 原子）平衡 200 ps 以获得该温度下的平衡晶格参数、平衡晶格夹角、平均原子坐标。
然后使用平衡结构参数建模扩胞，并在 NVT 系综下平衡 80 ps，然后在 NVE 系综下模拟 600 ps 采样。
从中抽取 300 帧提取二阶、三阶力常数。

以下必要的输入文件已包含在 `./TDEP/infiles_T300K` 目录中：

- infile.ucposcar：定义系统的原胞（晶格向量和平衡位置），格式为 VASP 的 POSCAR 格式。
- infile.ssposcar：定义系统的超胞（晶格向量和平衡位置），对应生成位置-力数据集的超胞，格式为 VASP 的 POSCAR 格式。
- infile.positions 和 infile.forces：分别包含位置数据和力数据。
- infile.meta：包含某些相关的元数据。
- infile.stat：包含分子动力学模拟输出的其它数据。

### 提取力常数

在包含上述所有输入文件的目录中，可以使用以下命令提取二阶有效力常数（IFCs）：

```bash
extract_forceconstants -rc2 100 -s 50 > fc2.log
```

此命令将通过最小二乘法拟合，从 infile.{positions,forces} 文件中找到最符合位置-力数据的二阶 IFCs 集合。
在拟合之前，独立的 IFC 数量将通过强制施加适当的对称性来减少。
这些对称性包括晶体结构的对称性，以及一般的平移和旋转不变性。
这通常会显著减少需要拟合的 IFC 数量。

几点说明：

- 选项 -rc2 X：指定二阶 IFC 的相互作用截断距离，即仅考虑距离小于 X 的原子对之间的相互作用。
- 指定较大的数值（例如 -rc2 100）将强制使用 infile.ssposcar 中适合的最大可能截断距离。指定较小的数值（例如 -rc2 0）将仅考虑最近邻之间的相互作用。
- 截断距离的重要性：这是一个需要检查的重要收敛参数。在完成本教程和其他相关教程后，建议重新计算不同 rc2 值下的声子色散关系等，体会截断距离的影响。
- 选项 -s 50：指定仅使用 infile.positions 和 infile.forces 中每 50 个样本中的一个。这控制了拟合中使用的样本数量。样本数量的重要性：样本数量也是一个重要的收敛参数。请注意，沿分子动力学轨迹密集采样的样本往往具有高度相关性。

请查看 extract_forceconstant 的输出结果，在本例中重定向到了文件 fc2.log。其中打印了许多有用的信息，尝试找到以下内容：

- 当前结构和超胞的最小（最近邻距离）和最大截断距离。
- 拟合的（不可约）二阶 IFC 的总数。
- 拟合中使用的总方程数（你可以思考如何计算这个数值？）。
- 拟合的 R² 值，这是衡量所获得的有效谐波 IFC 对位置-力数据描述效果的指标。

另外，运行结束后会生成两个输出文件：

- outfile.forceconstant：包含完整的二阶 IFC。
- outfile.irrifc_secondorder：包含不可约的二阶 IFC 列表。

### 声子色散关系与态密度 (DOS)

在进一步处理之前，需要将 outfile.forceconstant 复制或创建符号链接为 infile.forceconstant：

```bash
ln -s outfile.forceconstant infile.forceconstant
```

使用以下命令生成声子色散关系：

```bash
phonon_dispersion_relations
```

这将生成一系列输出文件，其中最重要的是 outfile.dispersion_relations，它包含色散关系的数据。
同时，还会生成一个 .gnuplot 文件，便于可视化声子色散关系。
运行以下命令绘制色散关系：

```bash
gnuplot --persist outfile.dispersion_relations.gnuplot
```

通过与文献比较，确认生成的声子色散关系是否与硅 (Si) 的预期声子色散关系一致。

如果需要生成声子态密度，可以添加 --dos 选项：

```
bash
phonon_dispersion_relations --dos
```

此命令同样会生成一个 gnuplot 脚本，用于查看 DOS：

```bash
gnuplot --persist outfile.phonon_dos.gnuplot
```

注意事项

- 布里渊区路径 (BZ path)：

如果未指定布里渊区路径（如上述例子），将使用当前空间群的默认路径。
若需自定义路径，可提供一个 infile.qpoints_dispersion 文件，并使用 --readpath 选项。
文件格式详见 TDEP 文件格式说明。
影响 DOS 的因素：

- 布里渊区网格、展宽方法及相关参数会影响生成的态密度。

这些参数可通过 phonon_dispersion_relations 的选项进行调整。
可运行以下命令查看所有选项，并花些时间研究适合的设置：

```bash
phonon_dispersion_relations -h
```

### 声子自由能

通过 phonon_dispersion_relations 命令的 --temperature（单一温度）或 --temperature_range（温度范围）选项，可以计算声子自由能。
例如，要获取在 0 到 1000 K 范围内的 50 个温度点的声子自由能，可以运行以下命令：

```bash
phonon_dispersion_relations --temperature_range 0 1000 50
```

此命令将生成一个 outfile.free_energy 文件，其中包含以下数据：温度、声子自由能、振动熵、热容。
要在 gnuplot 中绘制声子自由能图，例如运行以下命令：

```bash
plot "outfile.free_energy" u 1:2 w l
```

### 高阶力常数 (IFCs)

extract_forceconstants 也可以用来提取高阶 IFCs。
例如，要在截断距离为 4 Å 的条件下提取三阶 IFC，同时保留二阶 IFC，可以运行：

```bash
extract_forceconstants -rc2 100 -rc3 4 -s 50 > fc2_fc3.log
```

除了生成 outfile.forceconstant（二阶 IFC 文件），还会生成 outfile.forceconstant_thirdorder，其中包含三阶 IFC。

在日志文件中会打印一些有用的信息：
- 拟合了多少三阶 IFC？
- 通过在拟合中加入三阶 IFC，R² 是否有所提升？

### TDEP 的随机采样：sTDEP
简单来说，sTDEP 的思想是通过在（简谐）正则系综中使用力常数来近似原子位移分布，并利用系统中的真实力对该近似进行迭代改进，从而自洽地生成力常数。

与直接通过分子动力学（MD）模拟获取真实的核分布（代价更高）相比，sTDEP 方法通过采样近似的（有效）简谐分布，并在每次迭代后自洽更新力常数，直到收敛为止。

### 非谐声子计算

通过引入谱线形的影响，将其加入到声子色散关系中，从而可视化声子带的展宽和可能的混合。
可用于将模拟中得到的声子色散关系与非弹性中子散射实验进行对比。

要运行 --path 模式，可以执行以下命令：

```bash
mpirun /path/to/tdep/bin/lineshape --path -qg 3 3 3 --temperature 100 > path.log
```

命令参数：
- `--path`：指定使用路径计算模式。
- `-qg 3 3 3`：设置 q 点网格大小为 3×3×3（可根据需求后续进行收敛测试）。
- `--temperature 100`：设定计算温度为 100K。
- `> path.log`：将输出重定向到日志文件 path.log。

注意事项：
- 起始温度选择最低温度（如 100K）。
- 将 TDEP 的二进制文件路径替换为实际安装位置。
- 开始时可使用较小的 q 点网格，后续进行收敛性测试。

该模式允许沿布里渊区定义特定路径，默认路径与声子色散关系计算中的路径一致（即晶体高对称点之间的路径）。
通过此模式，可以清晰地展示声子谱线形对色散关系的影响，同时为实验数据对比提供基础。

通过使用 --readpath 选项，可以在计算中指定自定义的布里渊区路径，TDEP 将从 infile.qpoints_dispersion 文件中读取 q 点路径。
此外，可以通过 -nq（或 --nq_on_path）标志调整每两个高对称点之间的 q 点数量，默认值为 100，用于生成更密集的网格。

以下是 infile.qpoints_dispersion 的两个示例：

```
FCC                         ! Bravais lattice type
  100                       ! Number of points on each path
    4                       ! Number of paths between special points
GM  X                       ! Starting and ending special point
X   U                       !
K   GM                      !
GM  L                       !
```

```
CUSTOM                      !
  100                       ! Number of points on each path
    4                       ! Number of paths between special points
0.000 0.000 0.000   0.000 0.500 0.500 GM X
0.000 0.500 0.500   0.000 0.625 0.375 X  U
0.375 0.750 0.375   0.000 0.000 0.000 K  GM
0.000 0.000 0.000   0.000 0.500 0.000 GM L
```
