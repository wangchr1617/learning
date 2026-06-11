
# Packmol 入门

## 安装

将[官网](https://m3g.github.io/packmol/download.shtml) 上下载到的 `packmol.tar.gz` 上传到集群并解压。
进入解压后的目录，执行 `./configure` 配置编译文件，然后输入 `make` 编译即可。
编译完成后当前目录就出现了名为 `packmol` 的可执行文件。
最后修改 `~/.bashrc`：
```
export PATH=/home/changruiwang-ICME/Software/packmol-20.14.2:${PATH}
```
将其添加到环境变量即可。
输入 `packmol` 命令检查是否安装成功。

将 `.inp` 文件重定向到 `packmol` 即可运行程序，例如：
```
packmol < test.inp
```

---

## 简单示例

Packmol 程序需要预先提供输入文件和结构文件。

以生成液相 GeTe 为例，结构文件 `model.xyz` 如下所示：
```
2
Lattice="-2.92042768e-08 2.9909695976823745 2.99096938925063 2.990969671591604 -9.98001108e-08 2.990969459364865 2.9909696608638106 2.99096965747627 -8.91261311e-08" Properties=species:S:1:pos:R:3 pbc="T T T"
Ge      -0.00000000       0.00000000       0.00000000
Te       2.99096965       2.99096958       2.99096938
```

输入文件 run.inp 如下所示：
```
tolerance 2 # 要求原子间距不能小于 2.0 埃
output output.xyz # 指定输出文件
seed -1 # 使用随机种子生成随机结构
filetype xyz # 输出文件的格式
structure model.xyz # 指定输入文件
  number 150 # 分子数
  inside cube 0. 0. 0. 20.88330445 # [xmin ymin zmin d] 控制立方盒子的大小
end structure
```

将 model.xyz 和 run.inp 都放到当前目录，运行 `packmol < run.inp` 即可。

如果你设定的空间范围约束比较简单，Packmol 很快就会生成最终结构。
空间约束设定得越多、越复杂，生成结构需要的循环次数越多，总耗时就越高。

一般情况下，用 Packmol 构建的体系是比较松散的，其密度显著低于真实密度。
因此需要在 NPT 系综下跑一段分子动力学模拟调节盒子尺寸。