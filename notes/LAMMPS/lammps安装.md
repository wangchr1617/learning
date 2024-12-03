
# LAMMPS 安装

## LAMMPS make 安装

安装 LAMMPS 需要 MPICH、FFTW、以及 LAMMPS 安装包。
**MPICH** 在集群 Intel 编译器中有，可以通过 `which mpirun` 查看路径。
**FFTW** 在集群的 `/opt/software/` 目录下。
LAMMPS 安装包可以到 [LAMMPS 官网](https://note.youdao.com/) 下载。

1.  解压 LAMMPS 安装包：`tar -xzvf lammps-stable.tar.gz`；
2.  进入 LAMMPS 安装目录，并修改配置文件：

<!---->

    cd /home/zfli-ICME/software/lammps/lammps-29Sep2021/
    vim src/MAKE/OPTIONS/Makefile.g++_mpich                # 修改该文件
        MPI_INC =       -DMPICH_SKIP_MPICXX -DOMPI_SKIP_MPICXX=1
        MPI_PATH =  -L/opt/intel2020/intelpython3/bin/mpirun
        MPI_LIB =
        
        FFT_INC = -I/opt/software/fftw3.3.4/include
        FFT_PATH = -L/opt/software/fftw3.3.4/lib
        FFT_LIB =

3.  连接编译节点：`ssh cu02`；
4.  加载 Intel 编译环境：`module load intel/2020.1.217`；
5.  依次输入以下命令：

<!---->

    make yes-molecule
    make yes-kspace
    make yes-rigid
    make yes-manybody
    make g++_mpich

编译完成后 `src` 目录下会生成名为 `lmp\_g++\_mpich` 的可执行文件。

> 使用命令 `make package-installed` 查看编译了哪些包。

### LAMMPS-MLIP-interface 安装

安装 MTP 势函数的接口。

首先需要访问 [网址](https://gitlab.com/ashapeev/mlip-2) 下载 mlip 安装包。
我们仅编译需要用到的 mlip 接口库 **libinterface** 即可。
编译方法如下：

```
ssh cu02
module load intel/2020.1.217
cd mlip-2
mkdir build
cd build
cmake ..
make libinterface
```

编译完成后，该目录下会出现 **lib_mlip_interface.a** 文件，将其复制到 **interface-lammps-mlip-2-master** 目录下。

然后开始编译 MTP 势函数的接口。
先在 `preinstall.sh` 脚本末尾添加以下命令，在编译时将这些包也编译好：

```
make yes-molecule
make yes-kspace
make yes-rigid
make yes-manybody
make yes-REPLICA        # neb计算用
```

然后用下面的命令编译好接口文件，注意这里编译的是 **intel_cpu_intelmpi**，因为 **g++_mpich** 方式编译会报错：

```
module load intel/2020.1.217
./install.sh ../lammps-23Jun2022 intel_cpu_intelmpi
```

编译完成后该目录下会包含名为 **lmp_intel_cpu_intelmpi** 的可执行文件。

---

## LAMMPS cmake 安装

> cmake 和 make 安装方式不能混用，如果用 cmake 安装 LAMMPS，后续安装扩展包时不能用 make 安装。

### 检查并准备安装环境

首先检查 CMAKE，至少需要 3.10 及以上版本。

```
which cmake
cmake --version
```

cmake 方式安装可以分为三部分：

> 1. configuration：检测那些特性和设置应该被启用，以及如何编译 LAMMPS
> 2. compilation：生成和编译所有必要的源文件，构建库和可执行文件
> 3. installation：复制编译好的可执行文件到指定位置，不在依赖于源文件

配置和编译必须在单独的 build 目录中，并且该目录必须不同于源目录。

此外，源目录 `src` 必须保持原始状态，因此在此之前不能用make方式进行编译，因为make是在源目录中进行的。
可以用 `make no-all purge` 删除 `make` 创建的源文件。

最好是在顶层目录(lammps下)创建独立的目录，可以根据不同情况自己创建，比如编译并行版时为 build-parallel，编译串行版时为 build-serial。
这样在一次编译过程中所有的辅助文件都在该目录下。

### 安装

进入build目录，然后运行 `cmake ../cmake` 命令，屏幕上会输出配置进度，后面包括启用的特性，选项和编译器设置。

    tar -xzvf lammps-stable.tar.gz

    cd lammps-23Jun2022
    make build && cd build
    cmake ../cmake                      # 这一步并不会编译代码

    # 这两种方式任选一种都可以编译，一定是在build目录下，make方式还需要make install
    cmake --build .                     # 这一步才开始编译

    make -j 8
    make install

***

以下是运行 `cmake ../cmake` 命令后的输出，特别注意 **Could NOT find** 之后的内容，如果是必须的库，就需要先手动安装好。

### 设置选项

在编译时启动，关闭，调整某些设置可以用 **-D** 选项。格式为 `-D VARIABLE=VALUE` 。value一般为布尔值（on/off、yes/no、1/0）或字符串（字符串中若有空格需要用双引号括起来）

设置选项可以分类两类：

*   1\)cmake自带的设置选项
*   2\)lammps指定的设置选项

如果要编译 LAMMPS 的一些扩展包，设置方式为 `-D PKG_<NAME>=on` ；如果要禁用某个扩展包：`-D PKG_<NAME>=no` 。name即为扩展包的名称。

#### presets设置

当安装的扩展包很多的时候，都写在命令行中会很臃肿，且容易出错。可以通过设置preset fils，将需要启用的设置放到源目录的cmake/presets文件夹中，在命令行中指定该文件即可。preset文件通过 **-C** 加载。

    cmake -C ../cmake/presets/basic.cmake -D PKG_MISC=on ../cmake
    cmake -C ../cmake/presets/clang.cmake -C ../cmake/presets/most.cmake ../cmake
    cmake -C ../cmake/presets/basic.cmake -D BUILD_MPI=off ../cmake

#### 扩展包安装

##### VORONOI package

在DOWNLOAD\_VORO=yes时会到指定链接下载，但是不要联网，不过提示中会出现要下载的链接地址，手动下载后将压缩包复制到build目录下的 **build-common/voro\_build-prefix/src/voro++-0.4.6.tar.gz** 中，然后正常编译。

    cmake -C ../cmake/presets/basic.cmake -D PKG_VORONOI=yes -D DOWNLOAD_VORO=yes ../cmake
    cmake -C ../cmake/presets/basic.cmake -D PKG_VORONOI=yes -D DOWNLOAD_VORO=yes -D PKG_REPLICA=yes -D PKG_ML-SNAP=yes ../cmake

##### PYTHON package

要想将 LAMMPS 与 Python 连接起来，需要将 LAMMPS 通过共享库的方式编译（`BUILD_SHARED_LIBS`），而不是静态链接库。编译完之后会有 `liblammps.so` 动态链接库。

*   Python\_EXECUTABLE：Python 解释器位置
*   PKG\_PYTHON：是否安装 Python 包
*   BUILD\_SHARED\_LIBS：以共享连接库的方式安装

> 注：在和 `VORONOI` 扩展包一起编译时，`VORONOI` 默认会把生成的动态链接库删除，这时候需要进入 `voro_build-prefix/src/` 目录将压缩包解压，然后将 `voro_build-prefix/src/voro++-0.4.6/Makefile` 中最后删除动态链接库的命令注释掉即可。

    cmake -C ../cmake/presets/basic.cmake -D PKG_VORONOI=yes -D DOWNLOAD_VORO=yes -D PKG_REPLICA=yes -D PKG_ML-SNAP=yes -D PKG_PYTHON=yes -D Python_EXECUTABLE=/home/zfli-ICME/software/anaconda3/bin/python3.8 -D BUILD_SHARED_LIBS=yes ../cmake
    cmake --build .

    # 编译完成后需要再build目录下设置，有以下方式
    # 第一种，但是这种方式需要联网好像
    make install-python

    # 第二种，直接配置，in place use。PYTHONPATH指向根目录下的python目录；LD_LIBRARY_PATH指向build目录，是为了找到动态链接库
    # 采用这种方式的话，在编译时最好不要用`LAMMPS_MACHINE=name`选项，因为产生的`liblammps_mpi.so`动态链接库不会被解析，原始的python目录没有指定name
    # 将下面的两行添加到.bashrc文件，注意添加到conda初始化命令后面
    export PYTHONPATH=/home/zfli-ICME/software/lammps-23Jun2022/python:${PYTHONPATH}    # 或者 /home/zfli-ICME/software/lammps-23Jun2022/build-python-1/python/lib
    export LD_LIBRARY_PATH=/home/zfli-ICME/software/lammps-23Jun2022/build-python-1:${LD_LIBRARY_PATH}

---

## GPU版本 LAMMPS cmake 安装

    ssh gpu01

    # 需要自己将cmake路径添加到PATH环境变量中
    export PATH=/opt/software/cmake-3.15.1/bin/:$PATH

    mkdir build-gpu
    cmake -C ../cmake/presets/basic.cmake -D PKG_GPU=on -D GPU_API=cuda -D GPU_PREC=mixed -D GPU_ARCH=sm_86 -D CUDPP_OPT=yes ../cmake

    make -j 8
    make install            # 完成后会在该目录下生成一个lmp可执行文件，可以随意命名

    # 调用gpu版本的lammps做计算，需要加 -sf -pk 两个参数。-sf指在所有支持gpu加速的脚本命令前加上gpu前缀；-pk gpu后跟的是节点数，几块gpu就填几
    mpirun -np 8 lmp_gpu -sf gpu -pk gpu 1 -in in.file

<!---->

    cmake -c ../cmake/presets/kokkos-cuda.cmake -C ../cmake/presets/basic.cmake -D PKG_GPU=on -D GPU_API=cuda -D GPU_PREC=mixed -D GPU_ARCH=sm_86 -D CUDPP_OPT=yes -D PKG_ML-SNAP=yes ../cmake

