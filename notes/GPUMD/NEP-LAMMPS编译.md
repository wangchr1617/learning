
# NEP-LAMMPS编译

首先，`ssh cu02` 切换到编译节点。

然后加载 Intel 编译模块：

```
module load intel/2020.1.217
```

下载 Lammps 并解压，重命名为 lammps-nep：

```
git clone --depth=1 https://github.com/lammps/lammps
unzip lammps-release.zip
mv lammps-release lammps-nep
```

然后安装 NEP_CPU：
```
git clone https://github.com/brucefan1983/NEP_CPU.git
unzip NEP_CPU-main.zip
cd NEP_CPU-main/
cp src/* interface/lammps/USER-NEP/
cp interface/lammps/USER-NEP/ ../lammps-nep/src/
```

最后编译 LAMMPS：

```
cd ~/soft/lammps-nep/
mkdir build
cd build/
/home/yidongli-ICME/soft/cmake/bin/cmake ../cmake \
  -DCMAKE_BUILD_TYPE=Release \
  -DPKG_NEP=ON
make -j$(nproc) VERBOSE=1
make install
```

验证 LAMMPS 安装：`$HOME/soft/lammps-nep/build/lmp -help`

编译完成后使用 `logout` 退回登录节点。