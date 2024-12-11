
# ALKEMIE 集群机器学习势环境搭建

首先，`ssh gpu04` 切换到 gpu 队列的编译节点。

如果没有管理员权限的话，可以提前创建一个临时目录供编译使用：

```
mkdir -p $HOME/tmp
chmod 1777 $HOME/tmp
export TMPDIR=$HOME/tmp
```

## 前置安装

### 创建环境

考虑到集群不允许联网，使用命令 `conda create -n mlp --clone base` 通过克隆 base 环境创建虚拟环境 mlp。
然后使用命令 `conda activate mlp` 切换虚拟环境。

### 配置 GCC

集群 3 已经有编译好的 GCC，设置环境变量加载即可。

```
export PATH=/opt/gcc/gcc-9.3.0/bin:${PATH}
export LIBRARY_PATH=/opt/gcc/gcc-9.3.0/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/opt/gcc/gcc-9.3.0/lib64:$LD_LIBRARY_PATH
export C_INCLUDE_PATH=/opt/gcc/gcc-9.3.0/include:$C_INCLUDE_PATH
export CPLUS_INCLUDE_PATH=/opt/gcc/gcc-9.3.0/inclde:$CPLUS_INCLUDE_PATH
```

### 安装 GNU Make

集群安装的 GNU Make 版本是 3.82，这里给出升级到 4.4 版本的操作流程。

```
wget https://mirrors.tuna.tsinghua.edu.cn/gnu/make/make-4.4.tar.gz
tar -xvzf make-4.4.tar.gz
cd make-4.4
./configure --prefix=$HOME/soft/make
make -j$(nproc)
make install
export PATH=$HOME/soft/make/bin:$PATH
alias gmake='~/soft/make/bin/make'
export MAKE=$HOME/soft/make/bin/make
```

### 安装 Cmake

集群安装的 CMake 最高版本是 3.18.5，这里给出升级到 3.23.1 版本的操作流程。

```
wget https://github.com/Kitware/CMake/releases/download/v3.23.1/cmake-3.23.1-linux-x86_64.sh
sh ./cmake-3.23.1-linux-x86_64.sh --prefix=$HOME/soft/cmake/ --exclude-subdir
export PATH=$HOME/soft/cmake/bin:${PATH}
```

### 安装 GLIBC （似乎不是必要的）

使用 `ldd --version` 查看当前加载的 GLIBC 版本。
集群 3 安装的 GLIBC 版本是 2.17，可以手动升级到 2.34 版本。

```
wget https://ftp.gnu.org/gnu/glibc/glibc-2.34.tar.gz
tar -xvzf glibc-2.34.tar.gz
cd glibc-2.34
mkdir build
cd build
CFLAGS="-Wno-error" ../configure --prefix=$HOME/soft/glibc --with-ld=/opt/gcc/gcc-9.3.0/bin/ld --with-gcc=/opt/gcc/gcc-9.3.0/bin/gcc --with-g++=/opt/gcc/gcc-9.3.0/bin/g++
make -j1
make install
```

如果提示 `LD_LIBRARY_PATH shouldn't contain the current directory` 的话，
可以执行 `export LD_LIBRARY_PATH=$(echo $LD_LIBRARY_PATH | sed 's/:*\.*//g')` 避免之。

注意，不要全局设置 LD_LIBRARY_PATH 或 LD_PRELOAD，因为这会干扰系统的运行。
如果需要运行某些程序时使用自定义的 GLIBC，可以用以下命令临时设置：

```
LD_LIBRARY_PATH=$HOME/soft/glibc/lib:$LD_LIBRARY_PATH <command>
```

### 配置 Intel

集群 3 已经有编译好的 Intel 2020，设置环境变量加载即可。

加载 Intel 配置如下：

```
source /opt/intel2020/compilers_and_libraries_2020/linux/bin/compilervars.sh -arch intel64  
source /opt/intel2020/mkl/bin/mklvars.sh intel64
source /opt/intel2020/impi/2019.7.217/intel64/bin/mpivars.sh
```

或者指定环境变量如下所示：

```
export MANPATH=/opt/intel2020/compilers_and_libraries_2020.1.217/linux/mpi/man:$MANPATH
export MKLROOT=/opt/intel2020/compilers_and_libraries_2020.1.217/linux/mkl:$MKLROOT
export LD_LIBRARY_PATH=$MKLROOT/lib:$LD_LIBRARY_PATH
```

### 安装 cuDNN

访问 [网址](https://developer.nvidia.com/cudnn-downloads) 下载 cuDNN 安装包 cudnn-linux-x86_64-8.5.0.96_cuda11-archive.tar.xz。

```
tar -xvf cudnn-linux-x86_64-8.5.0.96_cuda11-archive.tar.xz
mv cudnn-linux-x86_64-8.5.0.96_cuda11-archive/ cuda

export CUDA_HOME=/usr/local/cuda-11.7
export PATH=$CUDA_HOME/bin:$PATH
export LD_LIBRARY_PATH=$CUDA_HOME/lib64:$LD_LIBRARY_PATH
export CUDNN_INCLUDE_PATH=~/soft/cuDNN/include
export CUDNN_LIBRARY_PATH=~/soft/cuDNN/lib
export LD_LIBRARY_PATH=$CUDNN_LIBRARY_PATH:$LD_LIBRARY_PATH

ls $CUDNN_INCLUDE_PATH/cudnn.h
ls $CUDNN_LIBRARY_PATH/libcudnn*
```

### 安装 NCCL

```
unzip nccl-2.18.zip
cd nccl-2.18/
make CUDA_HOME=/usr/local/cuda-11.7 PREFIX=$HOME/soft/nccl -j$(nproc)
make install PREFIX=$HOME/soft/nccl

export NCCL_HOME=$HOME/soft/nccl
export LD_LIBRARY_PATH=$NCCL_HOME/lib:$LD_LIBRARY_PATH
export CPATH=$NCCL_HOME/include:$CPATH
export PKG_CONFIG_PATH=$NCCL_HOME/lib/pkgconfig:$PKG_CONFIG_PATH
```

### 安装 Pytorch

访问 GitHub 下载 Pytorch 安装包 pytorch-v1.13.0.tar.gz 并解压。

```
conda install -c https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/main mkl mkl-include --no-deps --use-local
```

如果安装失败，访问 [网址](https://anaconda.org/anaconda/mkl/files?version=2022.0.1) 下载 mkl 和 mkl-include 库。
然后使用命令 `conda install mkl-2022.0.1-h06a4308_117.tar.bz2 --use-local` 安装。

```
untgz pytorch-v2.1.2.tar.gz
cd pytorch-v2.1.2/
mkdir -p build
cd build
cmake \
  -DBUILD_PYTHON=ON \
  -DBUILD_TEST=OFF \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX=$HOME/anaconda3/envs/mlp/lib/python3.9/site-packages \
  -DCMAKE_PREFIX_PATH=${CONDA_PREFIX:-$(dirname $(which conda))/../} \
  -DGLIBCXX_USE_CXX11_ABI=1 \
  -DCUDA_TOOLKIT_ROOT_DIR=$CUDA_HOME \
  -DCUDNN_INCLUDE_PATH=$CUDNN_INCLUDE_PATH \
  -DCUDNN_LIBRARY_PATH=$CUDNN_LIBRARY_PATH \
  -DNCCL_INCLUDE_DIR=$HOME/soft/nccl/include \
  -DNCCL_LIB_DIR=$HOME/soft/nccl/lib \
  -DUSE_NUMPY=ON \
  -DUSE_CUDA=ON \
  -DUSE_KINETO=OFF \
  -DBUILD_CAFFE2=OFF \
  -DUSE_MKL=ON \
  -DMKL_INCLUDE_DIR=$CONDA_PREFIX/include \
  -DMKL_LIBRARY=$CONDA_PREFIX/lib \
  -DPYTHON_EXECUTABLE=$HOME/anaconda3/envs/mlp/bin/python \
  ..
make -j$(nproc) VERBOSE=1
ls -l build/lib/libtorch_global_deps.so
```

如果库文件存在，说明构建成功。
构建完成后，返回到 PyTorch 的根目录，并使用 setup.py 进行开发模式安装。

```
export _GLIBCXX_USE_CXX11_ABI=1
export USE_CUDA=1
export TORCH_CUDA_ARCH_LIST="8.6"
export NCCL_INCLUDE=$HOME/soft/nccl/include
export NCCL_LIB=$HOME/soft/nccl/lib
export LD_LIBRARY_PATH=$(pwd)/build/lib:$NCCL_LIB:$LD_LIBRARY_PATH
export CXX=/opt/gcc/gcc-9.3.0/bin/g++
export CC=/opt/gcc/gcc-9.3.0/bin/gcc
export GIT_DISCOVERY_ACROSS_FILESYSTEM=1

python setup.py develop --verbose > build.log 2>&1
tail -n 20 build.log
```

如果安装成功，最后的日志会显示 `Installed /path/to/torch`.

最后，验证安装：

```
python
import torch
print(torch._C._GLIBCXX_USE_CXX11_ABI)
torch.cuda.nccl.is_available()
```

---

## 安装 Allegro、Deepmd-kit、NEP_CPU 及 Lammps

### 下载 lammps
git clone --depth=1 https://github.com/lammps/lammps
unzip lammps-release.zip

### 安装 Allegro
pip install ase-3.22.0-py3-none-any.whl
pip install opt_einsum-3.4.0-py3-none-any.whl
pip install opt_einsum_fx-0.1.4-py3-none-any.whl
pip install e3nn-0.5.4-py3-none-any.whl
pip install torch_runstats-0.2.0-py3-none-any.whl
pip install torch_ema-0.3-py3-none-any.whl

git clone --depth 1 https://github.com/mir-group/nequip
unzip nequip-0.6.1.zip
cd nequip-0.6.1
pip install .

git clone --depth 1 https://github.com/mir-group/allegro.git
unzip allegro-main.zip
cd allegro-main/
pip install .

git clone --depth 1 https://github.com/mir-group/pair_allegro.git
unzip pair_allegro-main.zip
cd pair_allegro-main/
bash patch_lammps.sh ../lammps-release/

---

### 安装 Deepmd-kit

git clone https://github.com/deepmodeling/deepmd-kit.git deepmd-kit
untgz deepmd-kit-3.0.0.tar.gz
cd deepmd-kit-3.0.0/
deepmd_source_dir=`pwd`

alias clientin='./ladder -u zb2201110 -p 19981019zhangdi auth'
alias clientout='./ladder -u zb2201110 auth --logout'

export DP_VARIANT=cuda
export DP_ENABLE_TENSORFLOW=0
export DP_ENABLE_PYTORCH=1

pip install -i https://pypi.tuna.tsinghua.edu.cn/simple . --timeout=1000 -v

HOROVOD_WITHOUT_GLOO=1 HOROVOD_WITH_PYTORCH=1 HOROVOD_GPU_OPERATIONS=NCCL HOROVOD_NCCL_HOME=/path/to/nccl pip install horovod mpi4py

cd $deepmd_source_dir/source
mkdir build
cd build
cmake -DENABLE_PYTORCH=TRUE -DCMAKE_PREFIX_PATH=/path/to/libtorch -DCMAKE_INSTALL_PREFIX=$deepmd_root -DUSE_CUDA_TOOLKIT=On ..
make -j$(nproc)
make install
ls $deepmd_root/bin
ls $deepmd_root/lib

### 安装 NEP_CPU

```
git clone https://github.com/brucefan1983/NEP_CPU.git
cp NEP_CPU/src/* NEP_CPU/interface/lammps/USER-NEP/
cp NEP_CPU/interface/lammps/USER-NEP/ YOUR_LAMMPS_PATH/src/
```

### 安装 Lammps

cd lammps
rm -rf build
mkdir build
cd build
cmake ../cmake \
  -D CMAKE_BUILD_TYPE=Release \
  -D CMAKE_PREFIX_PATH=/content/libtorch \
  -D MKL_INCLUDE_DIR=$MKLROOT/include \
  -D MKL_LIBRARIES=$MKLROOT/lib/intel64 \
  -D PKG_PLUGIN=ON \
  -D PKG_NEP=on \
  -D LAMMPS_INSTALL_RPATH=ON \
  -D BUILD_SHARED_LIBS=yes \
  -D CMAKE_INSTALL_PREFIX= \ # 指定安装路径
  -D CMAKE_INSTALL_LIBDIR=lib
make -j$(nproc)
make install