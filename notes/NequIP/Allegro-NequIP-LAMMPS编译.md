
# Allegro - NequIP - LAMMPS 编译

首先，`ssh gpu04` 切换到 gpu 队列的编译节点。

### MLP 环境搭建

考虑到集群不允许联网，使用命令 `conda create -n mlp --clone base` 通过克隆 base 环境创建虚拟环境 mlp。
然后使用命令 `conda activate mlp` 切换虚拟环境。

> 如果需要删除环境，使用 `conda remove --name mlp --all` 删除虚拟环境。

### CUDA 环境

设置 CUDA 相关的环境变量如下：

```
export PATH=/usr/local/cuda-11.7/bin:$PATH
export LD_LIBRARY_PATH=/usr/local/cuda-11.7/lib64:$LD_LIBRARY_PATH
export CUDA_HOME=/usr/local/cuda-11.7
export CUDA_ROOT=/usr/local/cuda-11.7
```

### 安装 Pytorch

访问 [网址](https://download.pytorch.org/whl/torch/) 下载 Pytorch 安装包 `torch-2.0.1+cu117.with.pypi.cudnn-cp39-cp39-linux_x86_64.whl`。

依次安装依赖如下所示：

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
pip install triton-2.0.0-1-cp39-cp39-manylinux2014_x86_64.manylinux_2_17_x86_64.whl --no-deps
pip install nvidia_cusparse_cu11-11.7.4.91-py3-none-manylinux1_x86_64.whl --no-deps
pip install nvidia_cuda_cupti_cu11-11.7.101-py3-none-manylinux1_x86_64.whl --no-deps
pip install nvidia_nccl_cu11-2.14.3-py3-none-manylinux1_x86_64.whl --no-deps
pip install nvidia_nvtx_cu11-11.7.91-py3-none-manylinux1_x86_64.whl --no-deps
```

然后安装 Pytorch 和 MKL INCLUDE：

```
pip install torch-2.0.1+cu117.with.pypi.cudnn-cp39-cp39-linux_x86_64.whl
pip install mkl_include-2024.2.2-py2.py3-none-manylinux1_x86_64.whl
```

最后，使用 `python3 -c "import torch; print(torch.__version__)"` 验证安装。

### 安装 libtorch

访问 [网址](https://download.pytorch.org/libtorch/cu117/) 下载 `libtorch-shared-with-deps-2.0.1%2Bcu117.zip`。

```
unzip libtorch-shared-with-deps-2.0.1%2Bcu117.zip
```

之所以不选择 cxx11-abi 版本是因为 `strings /usr/lib64/libstdc++.so.6 | grep GLIBCXX` 发现集群系统的 libstdc++.so.6 版本最高支持到 GLIBCXX_3.4.19，
而支持 cxx11-abi 的 PyTorch 通常需要至少 GLIBCXX_3.4.20 或更高版本。

### 安装 Allegro

依次安装依赖如下所示：

```
pip install ase-3.22.0-py3-none-any.whl
pip install opt_einsum-3.4.0-py3-none-any.whl
pip install opt_einsum_fx-0.1.4-py3-none-any.whl
pip install e3nn-0.5.4-py3-none-any.whl
pip install torch_runstats-0.2.0-py3-none-any.whl
pip install torch_ema-0.3-py3-none-any.whl
```

安装 NequIP：

```
git clone --depth 1 https://github.com/mir-group/nequip
unzip nequip-0.6.1.zip
cd nequip-0.6.1
pip install .
```

安装 Allegro：

```
git clone --depth 1 https://github.com/mir-group/allegro.git
unzip allegro-main.zip
cd allegro-main/
pip install .
```

### 安装 Openmpi

访问 [网址](https://www.open-mpi.org/software/ompi/v4.0/) 下载 Openmpi 的离线安装包，然后解压编译：

```
tar -xzvf openmpi-4.0.5.tar.gz
cd openmpi-4.0.5
./configure --prefix=$HOME/soft/openmpi
make
make install
```

设置 Openmpi 环境变量如下:

```
MPI_HOME=$HOME/soft/openmpi
export PATH=${MPI_HOME}/bin:$PATH
export LD_LIBRARY_PATH=${MPI_HOME}/lib:$LD_LIBRARY_PATH
export MANPATH=${MPI_HOME}/share/man:$MANPATH
```

### 安装 LAMMPS

首先下载 LAMMPS 稳定版的离线安装包并解压：

```
git clone --depth=1 https://github.com/lammps/lammps
unzip lammps-release.zip
mv lammps-release lammps-allegro
```

然后安装 pair_allegro-main：

```
git clone --depth 1 https://github.com/mir-group/pair_allegro.git
unzip pair_allegro-main.zip
cd pair_allegro-main/
./patch_lammps.sh ../lammps-allegro/
```

最后，编译安装 LAMMPS：

```
cd ~/soft/lammps-allegro/
mkdir build
cd build/
/home/yidongli-ICME/soft/cmake/bin/cmake ../cmake \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_CXX_COMPILER=$HOME/soft/openmpi/bin/mpicxx \
  -DCMAKE_CXX_FLAGS="-L$HOME/soft/openmpi/lib -lmpi" \
  -DCMAKE_C_COMPILER=$HOME/soft/openmpi/bin/mpicc \
  -DCMAKE_PREFIX_PATH=$HOME/soft/libtorch \
  -DMKL_INCLUDE_DIR=$HOME/anaconda3/envs/mlp/include \
  -DMPI_CXX_INCLUDE_PATH=$HOME/soft/openmpi/include \
  -DMPI_CXX_LIBRARIES=$HOME/soft/openmpi/lib/libmpi.so
make -j$(nproc) VERBOSE=1
make install
```

验证 LAMMPS 安装：`$HOME/soft/lammps-allegro/build/lmp -help`

编译完成后使用 `logout` 退回登录节点。