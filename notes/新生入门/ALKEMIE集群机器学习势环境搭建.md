
# ALKEMIE 集群机器学习势环境搭建

首先，`ssh gpu04` 切换到 gpu 队列的编译节点。

如果没有管理员权限的话，可以提前创建一个临时目录供编译使用：

```
mkdir -p $HOME/tmp
chmod 1777 $HOME/tmp
export TMPDIR=$HOME/tmp
```

## 前置安装

### 配置 GCC

集群 3 已经有编译好的 GCC 9.3，设置环境变量加载即可。

```
export PATH=/opt/gcc/gcc-9.3.0/bin:${PATH}
export LIBRARY_PATH=/opt/gcc/gcc-9.3.0/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/opt/gcc/gcc-9.3.0/lib64:$LD_LIBRARY_PATH
export C_INCLUDE_PATH=/opt/gcc/gcc-9.3.0/include:$C_INCLUDE_PATH
export CPLUS_INCLUDE_PATH=/opt/gcc/gcc-9.3.0/inclde:$CPLUS_INCLUDE_PATH
```

使用 `gcc --version` 确保 GCC 环境加载正确。

### 安装 GLIBC

使用 `ldd --version` 查看当前加载的 GLIBC 版本。
集群 3 安装的 GLIBC 版本是 2.17，可以手动升级到 2.34 版本。

```
wget https://ftp.gnu.org/gnu/glibc/glibc-2.34.tar.gz
tar -xvzf glibc-2.34.tar.gz
cd glibc-2.34/
mkdir build
cd build
../configure --prefix=$HOME/soft/glibc CFLAGS="-Og -g -g3 -ggdb -gdwarf-4" CXXFLAGS="-Og -g -g3 -ggdb -gdwarf-4" --disable-werror
make -j1 VERBOSE=1
make install
```

注意，`-Og` 是必须得指定的，否则会报 "glibc cannot be compiled without optimization" 错误；`--disable-werror` 也必须设置，否则会把警告视为错误，无法编译。
如果 `configure` 提示 `LD_LIBRARY_PATH shouldn't contain the current directory` 的话，
可以执行 `export LD_LIBRARY_PATH=$(echo $LD_LIBRARY_PATH | sed 's/:*\.*//g')` 避免之。

注意，不要全局设置 `LD_LIBRARY_PATH` 或 `LD_PRELOAD`，因为这会干扰系统的运行。
如果需要运行某些程序时使用自定义的 GLIBC，可以用以下命令临时设置：

```
LD_LIBRARY_PATH=$HOME/soft/glibc/lib:$LD_LIBRARY_PATH <command>
```

使用 `$HOME/soft/glibc/lib/ld-linux-x86-64.so.2 --version` 验证安装是否成功。
然后设置 glibc 环境变量，使其对所有编译工具生效：

```
export LD_LIBRARY_PATH=$HOME/soft/glibc/lib:$LD_LIBRARY_PATH
export LIBRARY_PATH=$HOME/soft/glibc/lib:$LIBRARY_PATH
```

### 安装 Cmake

集群安装的 CMake 最高版本是 3.18.5，如果需要升级版本，以下是升级到 3.23.1 版本的操作流程。

```
wget https://github.com/Kitware/CMake/releases/download/v3.23.1/cmake-3.23.1-linux-x86_64.sh
chmod 777 cmake-3.25.0-linux-x86_64.sh
mkdir -p $HOME/soft/cmake/
$HOME/soft/glibc/lib/ld-linux-x86-64.so.2 \
--library-path $HOME/soft/glibc/lib:/usr/lib64:/opt/gcc/gcc-9.3.0/lib64 \
/bin/bash ./cmake-3.25.0-linux-x86_64.sh --prefix=$HOME/soft/cmake/ --exclude-subdir

```

设置路径：`export PATH=$HOME/soft/cmake/bin:$PATH`

验证 CMake 安装：`$HOME/soft/cmake/bin/cmake --version`

### 安装 GNU Make

集群安装的 GNU Make 版本是 3.82，如果需要升级版本，以下是升级到 4.4 版本的操作流程。

```
wget https://mirrors.tuna.tsinghua.edu.cn/gnu/make/make-4.4.tar.gz
tar -xvzf make-4.4.tar.gz
cd make-4.4
$HOME/soft/glibc/lib/ld-linux-x86-64.so.2 \
--library-path $HOME/soft/glibc/lib:/opt/gcc/gcc-9.3.0/lib64:/usr/lib64:/lib64 \
/bin/bash ./configure --prefix=$HOME/soft/make
make -j$(nproc)
make install
```

设置路径：
```
export PATH=$HOME/soft/make/bin:$PATH
alias gmake='$HOME/soft/make/bin/make'
```

验证 make 安装：`gmake --version`

### 安装 Openmpi

访问 [网址](https://www.open-mpi.org/software/ompi/v4.1/) 下载 Openmpi 安装包并拖入集群。

```
export CUDA_HOME=/usr/local/cuda-11.7
export PATH=$CUDA_HOME/bin:$PATH
export LD_LIBRARY_PATH=$CUDA_HOME/lib64:$LD_LIBRARY_PATH

untgz cd openmpi-4.1.6.tar.gz
cd openmpi-4.1.6/
$HOME/soft/glibc/lib/ld-linux-x86-64.so.2 \
--library-path $HOME/soft/glibc/lib:/opt/gcc/gcc-9.3.0/lib64:$CUDA_HOME/lib64:/usr/lib64:/lib64 \
/bin/bash ./configure --prefix=$HOME/openmpi \
            CC=/opt/gcc/gcc-9.3.0/bin/gcc \
            CXX=/opt/gcc/gcc-9.3.0/bin/g++ \
            CFLAGS="-I$HOME/soft/glibc/include" \
            CXXFLAGS="-I$HOME/soft/glibc/include" \
            LDFLAGS="-Wl,-rpath,$HOME/soft/glibc/lib -L$HOME/soft/glibc/lib" \
            --with-cuda=$CUDA_HOME \
            --enable-mpi-cxx
$HOME/soft/make/bin/make -j$(nproc)
$HOME/soft/make/bin/make install
```

设置环境变量：
```
MPI_HOME=$HOME/openmpi
export PATH=${MPI_HOME}/bin:$PATH
export LD_LIBRARY_PATH=${MPI_HOME}/lib:$LD_LIBRARY_PATH
export MANPATH=${MPI_HOME}/share/man:$MANPATH
export OMPI_MCA_orte_tmpdir_base=$HOME/tmp
export OMPI_MCA_shmem_mmap_enable_nfs_warning=0
```

验证 OpenMPI 安装：`mpirun --version`

---

## 环境搭建

考虑到集群不允许联网，使用命令 `conda create -n mlp --clone base` 通过克隆 base 环境创建虚拟环境 mlp。
然后使用命令 `conda activate mlp` 切换虚拟环境。
如果需要删除环境，使用 `conda remove --name mlp --all` 删除虚拟环境。

### 安装 Pytorch

访问 [网址](https://download.pytorch.org/whl/torch/) 下载 Pytorch 安装包 `torch-2.0.1+cu117.with.pypi.cudnn-cp39-cp39-linux_x86_64.whl`。

```
pip install nvidia_cusolver_cu11-11.4.0.1-py3-none-manylinux1_x86_64.whl --no-deps
pip install nvidia_cuda_nvrtc_cu11-11.7.99-py3-none-manylinux1_x86_64.whl --no-deps
pip install nvidia_cublas_cu11-11.10.3.66-py3-none-manylinux1_x86_64.whl --no-deps
pip install nvidia_cuda_runtime_cu11-11.7.99-py3-none-manylinux1_x86_64.whl --no-deps
pip install nvidia_curand_cu11-10.2.10.91-py3-none-manylinux1_x86_64.whl --no-deps
pip install nvidia_cufft_cu11-10.9.0.58-py3-none-manylinux1_x86_64.whl --no-deps
pip install cmake-3.25.0-py2.py3-none-manylinux_2_12_x86_64.manylinux2010_x86_64.whl --no-deps
untgz lit-15.0.7.tar.gz
cd lit-15.0.7/
python setup.py install
pip install nvidia_cudnn_cu11-8.5.0.96-py3-none-manylinux1_x86_64.whl --no-deps
pip install ./triton-2.0.0-1-cp39-cp39-manylinux2014_x86_64.manylinux_2_17_x86_64.whl --no-deps
pip install nvidia_cusparse_cu11-11.7.4.91-py3-none-manylinux1_x86_64.whl --no-deps
pip install nvidia_cuda_cupti_cu11-11.7.101-py3-none-manylinux1_x86_64.whl --no-deps
pip install nvidia_nccl_cu11-2.14.3-py3-none-manylinux1_x86_64.whl --no-deps
pip install nvidia_nvtx_cu11-11.7.91-py3-none-manylinux1_x86_64.whl --no-deps
pip install torch-2.0.1+cu117.with.pypi.cudnn-cp39-cp39-linux_x86_64.whl
pip install mkl_include-2024.2.2-py2.py3-none-manylinux1_x86_64.whl
```

最后，使用 `python3 -c "import torch; print(torch.__version__)"` 验证安装。

### 安装 libtorch

访问 [网址](https://download.pytorch.org/libtorch/cu117/) 下载 `libtorch-cxx11-abi-shared-with-deps-2.0.1%2Bcu117.zip`。

```
unzip libtorch-cxx11-abi-shared-with-deps-2.0.1%2Bcu117.zip
```

---

## 安装 Allegro、Deepmd-kit、NEP_CPU 及 Lammps

### 下载 lammps

git clone --depth=1 https://github.com/lammps/lammps
unzip lammps-release.zip

### 安装 Allegro

```
pip install ase-3.22.0-py3-none-any.whl
pip install opt_einsum-3.4.0-py3-none-any.whl
pip install opt_einsum_fx-0.1.4-py3-none-any.whl
pip install e3nn-0.5.4-py3-none-any.whl
pip install torch_runstats-0.2.0-py3-none-any.whl
pip install torch_ema-0.3-py3-none-any.whl
```

```
git clone --depth 1 https://github.com/mir-group/nequip
unzip nequip-0.6.1.zip
cd nequip-0.6.1
pip install .
```

```
git clone --depth 1 https://github.com/mir-group/allegro.git
unzip allegro-main.zip
cd allegro-main/
pip install .
```

```
git clone --depth 1 https://github.com/mir-group/pair_allegro.git
unzip pair_allegro-main.zip
cd pair_allegro-main/
./patch_lammps.sh ../lammps-release/
```

### 安装 NEP_CPU

```
git clone https://github.com/brucefan1983/NEP_CPU.git
unzip NEP_CPU-main.zip
cd NEP_CPU-main/
cp src/* interface/lammps/USER-NEP/
cp interface/lammps/USER-NEP/ ../lammps-release/src/
```

### 安装 Lammps

```
cd ~/lammps-release/
mkdir build
cd build/
$HOME/soft/glibc/lib/ld-linux-x86-64.so.2 \
--library-path $HOME/soft/glibc/lib:/opt/gcc/gcc-9.3.0/lib64:/usr/local/cuda-11.7/lib64:/usr/lib64:/lib64 \
$HOME/soft/cmake/bin/cmake ../cmake \
  -D CMAKE_BUILD_TYPE=Release \
  -D CMAKE_PREFIX_PATH=$HOME/soft/libtorch \
  -D CMAKE_EXE_LINKER_FLAGS="-Wl,-rpath,$HOME/soft/glibc/lib" \
  -D CMAKE_SHARED_LINKER_FLAGS="-L$HOME/soft/glibc/lib" \
  -D CMAKE_CXX_STANDARD_LIBRARIES="$HOME/soft/glibc/lib/libc.so" \
  -D CUDA_CUDART_LIBRARY=/usr/local/cuda-11.7/lib64/libcudart.so \
  -D CUDA_TOOLKIT_ROOT_DIR=/usr/local/cuda-11.7 \
  -D MKL_INCLUDE_DIR=$HOME/anaconda3/envs/mlp/include \
  -D PKG_PLUGIN=ON \
  -D PKG_NEP=ON \
  -D LAMMPS_INSTALL_RPATH=ON \
  -D CMAKE_INSTALL_PREFIX=$HOME/soft/lammps-install \
  -D CMAKE_INSTALL_LIBDIR=lib \
  -D PKG_PYTHON=yes \
  -D Python_EXECUTABLE=$HOME/anaconda3/envs/mlp/bin/python3.9 \
  -D BUILD_SHARED_LIBS=yes \
  -D CMAKE_MAKE_PROGRAM=$HOME/soft/make/bin/make
$HOME/soft/make/bin/make -j$(nproc) VERBOSE=1
$HOME/soft/make/bin/make install
```

验证 LAMMPS 安装：`$HOME/soft/lammps-install/bin/lmp -help`

