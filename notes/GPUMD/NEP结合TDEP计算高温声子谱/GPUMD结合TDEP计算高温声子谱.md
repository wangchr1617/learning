
# GPUMD 结合 TDEP 计算高温声子谱

本教程将介绍如何使用 TDEP 从 GPUMD 模拟轨迹中提取有效的原子间力常数（IFCs），并进一步计算声子色散关系、态密度和振动自由能等。

### GPUMD 模拟平衡结构参数

为了准备 TDEP 的输入文件，我们需要通过分子动力学模拟得到不同温度下的平衡结构参数。
以 GeTe 菱方相原胞的 10 × 10 × 10 超胞（共 2000 个原子）作为初始的模型文件 model.xyz。

使用 Python 生成 GPUMD 的输入文件 run.in：

```
temp_list = [
    100, 200, 300, 400, 500, 600, 700, 800
]

content = f"""
potential   /home/changruiwang-ICME/GAP_NEP/train/train/nep.txt
velocity    {temp_list[0]}

time_step   1
"""

for temp in temp_list:
    content += f"""
ensemble    nvt_bdp {temp} {temp} 100
run         10000

ensemble    npt_scr {temp} {temp} 100 0 0 0 0 0 0 100 100 50 50 50 50 1000
run         20000

ensemble    npt_scr {temp} {temp} 100 0 0 0 0 0 0 100 100 50 50 50 50 1000
dump_thermo 100
dump_exyz   100 0 0
run         100000
"""

with open("run.in", "w") as file:
    file.write(content)
```

GPUMD 作业结束后，会输出 dump.xyz 轨迹文件。

使用 Python 提取轨迹文件：

```
import os

input_file = 'dump.xyz'
temp_list = [
    100, 200, 300, 400, 500, 600, 700, 800
]

natoms = 2000
frames_per_file = 1000
lines_per_file = frames_per_file * (natoms + 2)
output_folder = 'split_xyz_files'

if not os.path.exists(output_folder):
    os.makedirs(output_folder)

with open(input_file, 'r') as file:
    lines = file.readlines()

total_lines = len(lines)
num_files = total_lines // lines_per_file

for i in range(num_files):
    start_index = i * lines_per_file
    end_index = (i + 1) * lines_per_file
    output_file = os.path.join(output_folder, f'split_{temp_list[i]}.xyz')
    with open(output_file, 'w') as split_file:
        split_file.writelines(lines[start_index:end_index])

if total_lines % lines_per_file != 0:
    output_file = os.path.join(output_folder, f'split_{num_files+1}.xyz')
    with open(output_file, 'w') as split_file:
        split_file.writelines(lines[num_files*lines_per_file:])
```

然后使用 Python 从不同温度对应的轨迹文件中提取平衡结构参数并写入 .vasp 文件中：

```
from ase import Atoms
from ase.io import read
from ase.io.vasp import write_vasp
import numpy as np
import os

temp_list = [
    100, 200, 300, 400, 500, 600, 700, 800
]
output_folder = 'split_vasp_files'

if not os.path.exists(output_folder):
    os.makedirs(output_folder)

def wrap_to_unit_cell(atom_pos, cell_matrix):
    inv_cell_matrix = np.linalg.inv(cell_matrix)
    fractional_coords = np.dot(atom_pos, inv_cell_matrix)
    fractional_coords = fractional_coords % 1.0
    return fractional_coords

def calculate_order_parameter(atoms, prim_cell_matrix, central_element='Ge', neighboring_element='Te'):
    atoms.center()
    positions = atoms.get_positions()
    central_atoms = [i for i, element in enumerate(atoms.get_chemical_symbols()) if element == central_element]
    neighboring_atoms = [i for i, element in enumerate(atoms.get_chemical_symbols()) if element == neighboring_element]
    central_positions = wrap_to_unit_cell(positions[central_atoms], prim_cell_matrix)
    neighboring_positions = wrap_to_unit_cell(positions[neighboring_atoms], prim_cell_matrix)
    tau = (np.mean((neighboring_positions)) - np.mean((central_positions)) - 0.5 * np.ones(3)) % 1.0
    return tau

for temp in temp_list:
    frames = read(f"./split_xyz_files/split_{temp}.xyz", index=":")
    a, b, c, alpha, beta, gamma = [], [], [], [], [], []

    tau_list = []
    for atoms in frames:
        cellpar = atoms.cell.cellpar()
        a.append(cellpar[0] / 10)
        b.append(cellpar[1] / 10)
        c.append(cellpar[2] / 10)
        alpha.append(cellpar[3])
        beta.append(cellpar[4])
        gamma.append(cellpar[5])
        prim_cell_matrix = atoms.cell / 10
        tau = calculate_order_parameter(atoms, prim_cell_matrix, central_element='Ge', neighboring_element='Te')
        tau_list.append(tau)
    
    cellpar_avg = (np.mean(a), np.mean(b), np.mean(c), np.mean(alpha), np.mean(beta), np.mean(gamma))
    a = np.mean(cellpar_avg[:3])
    alpha = np.mean(cellpar_avg[3:])
    tau_average = np.mean(tau_list, axis=0)
    scaled_positions = np.mean(tau_average)
    
    atoms = Atoms('GeTe')
    atoms.set_cell((a, a, a, alpha, alpha, alpha), scale_atoms=True)
    
    positions = [
        [0.0, 0.0, 0.0],
        [scaled_positions, scaled_positions, scaled_positions]
    ]
    atoms.set_scaled_positions(positions)
    write_vasp(f"./{output_folder}/split_{temp}.vasp", atoms, direct=True, sort=True)
```

### GPUMD 模拟采样

然后平衡结构的 6 × 6 × 6 超胞模拟采样，首先在 NVT 系综下平衡 80 ps，然后在 NVE 系综下模拟 600 ps 采样 300 帧。
以 300 K 为例，使用 Python 生成 GPUMD 的输入文件 run.in：

```
temp_list = [
    300, 400, 500, 510, 520, 530, 540, 550,
    560, 570, 580, 590, 600, 650, 700, 800
]

for temp in temp_list:
	content = f"""
potential   /home/changruiwang-ICME/GAP_NEP/train/train/nep.txt
velocity    {temp}

time_step   1

ensemble    nvt_bdp {temp} {temp} 100
run         80000

ensemble    nve
dump_thermo 2000
dump_exyz   2000 0 1
run         600000
"""

	with open(f"./split_{temp}/run.in", "w") as file:
		file.write(content)
```

然后使用 Python 生成 TDEP 的输入文件：

```
import numpy as np

dump_file = 'dump.xyz'
thermo_file = 'thermo.out'

natoms = 432
lines_per_frame = natoms + 2
time_step_interval = 2000
custom_temperature = 300

with open(dump_file, 'r') as file:
    dump_lines = file.readlines()

with open(thermo_file, 'r') as file:
    thermo_lines = file.readlines()

total_frames_dump = len(dump_lines) // lines_per_frame
total_frames_thermo = len(thermo_lines)
if total_frames_dump != total_frames_thermo:
    raise ValueError("dump.xyz 和 thermo.out 的帧数不一致！")

with open('infile.positions', 'w') as pos_file, \
     open('infile.forces', 'w') as force_file, \
     open('infile.stat', 'w') as stat_file, \
     open('infile.meta', 'w') as meta_file:

    meta_file.write(f"{natoms}\n")
    meta_file.write(f"{total_frames_dump}\n")
    meta_file.write(f"{time_step_interval}\n")
    meta_file.write(f"{custom_temperature}\n")

    for frame in range(total_frames_dump):
        dump_start = frame * lines_per_frame
        header = dump_lines[dump_start + 1]
        atom_lines = dump_lines[dump_start + 2 : dump_start + 2 + natoms]

        lattice_str = header.split('Lattice="')[1].split('"')[0]
        lattice = np.array(list(map(float, lattice_str.split()))).reshape(3, 3)

        stress_str = header.split('stress="')[1].split('"')[0]
        stress = np.array(list(map(float, stress_str.split()))).reshape(3, 3)
        stress_components = [stress[0, 0], stress[1, 1], stress[2, 2], 
                             stress[2, 0], stress[2, 1], stress[0, 1]] 

        positions = []
        forces = []
        for atom_line in atom_lines:
            parts = atom_line.split()
            pos = np.array(list(map(float, parts[1:4])))
            force = list(map(float, parts[4:7])) 
            frac_pos = np.linalg.solve(lattice.T, pos) 
            positions.append(frac_pos)
            forces.append(force)

        pos_file.writelines(' '.join(map(str, frac_pos)) + '\n' for frac_pos in positions)
        force_file.writelines(' '.join(map(str, force)) + '\n' for force in forces)

        thermo_values = list(map(float, thermo_lines[frame].split()))
        T, K, U, Px, Py, Pz = thermo_values[:6]
        E = K + U 
        P = (Px + Py + Pz) / 3 
        time = (frame + 1) * time_step_interval 

        stat_file.write(
            f"{frame + 1} {time:.6f} {E:.6f} {U:.6f} {K:.6f} {T:.6f} {P:.6f} "
            + ' '.join(map(str, stress_components)) + '\n'
        )
```

然后重命名 POSCAR 和 SC666.vasp 为 infile.ucposcar 和 infile.ssposcar 即可。

检查以下必要的输入文件是否已包含在当前目录中：

- `infile.ucposcar`：定义系统的原胞（晶格向量和平衡位置），格式为 VASP 的 POSCAR 格式。
- `infile.ssposcar`：定义系统的超胞（晶格向量和平衡位置），对应生成位置-力数据集的超胞，格式为 VASP 的 POSCAR 格式。
- `infile.positions`：位置数据。
- `infile.forces`：力数据。
- `infile.meta`：包含某些相关的元数据。
- `infile.stat`：包含分子动力学模拟输出的其它数据。

### 提取力常数

在包含上述所有输入文件的目录中，可以使用以下命令提取原子间力常数（IFCs）：

```bash
extract_forceconstants -rc2 10.5 -rc3 8 -rc4 5
```

此命令将通过最小二乘法拟合，从 infile.{positions,forces} 文件中找到最符合位置-力数据的 IFCs 集合。
在拟合之前，独立的 IFC 数量将通过强制施加适当的对称性来减少。
这些对称性包括晶体结构的对称性，以及一般的平移和旋转不变性。
这通常会显著减少需要拟合的 IFC 数量。

选项 -rc2/3/4 X：指定 IFC 的相互作用截断距离，即仅考虑距离小于 X 的原子对之间的相互作用。截断距离是一个需要检查的重要收敛参数。
样本数量也是一个重要的收敛参数。请注意，沿分子动力学轨迹密集采样的样本往往具有高度相关性。

### 声子色散关系与态密度 (DOS)

在进一步处理之前，需要将 outfile.forceconstant 复制或创建符号链接为 infile.forceconstant：

```bash
cp outfile.forceconstant infile.forceconstant
cp outfile.forceconstant_thirdorder infile.forceconstant_thirdorder
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

### 非谐声子计算

通过引入谱线形的影响，将其加入到声子色散关系中，从而可视化声子带的展宽和可能的混合。
可用于将模拟中得到的声子色散关系与非弹性中子散射实验进行对比。

要运行 --path 模式，可以执行以下命令：

```bash
mpirun /path/to/tdep/bin/lineshape --path -qg 6 6 6 -ne 600 --temperature 300 > path.log
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

使用 Python 脚本可以方便地可视化非谐声子谱。

```
from matplotlib.colors import LogNorm
import h5py as h5
import numpy as np
import matplotlib.pyplot as plt

f = h5.File("./outfile.phonon_spectral_function.hdf5", "r")
print(f.keys())

x = np.array(f.get("q_values"))
y = np.array(f.get("energy_values"))
gz = np.array(f.get("spectral_function"))
xt = np.array(f.get('q_ticks'))
xl = f.attrs.get('q_tick_labels').split()
gz=gz+1E-2
gx, gy = np.meshgrid(x,y)

plt.pcolormesh(gx, gy, gz, norm=LogNorm(vmin=gz.min(), vmax=gz.max()), cmap='afmhot')
plt.axis([x.min(), x.max(), y.min(), 5])
plt.xlabel("KPATH")
plt.ylabel("Frequency (THz)")
plt.xticks(xt, xl)

plt.savefig("phon.png", bbox_inches='tight')
```
