
# ALKEMIE 集群 3 机器学习势环境搭建

首先，`ssh gpu04` 切换到 gpu 队列的编译节点。
编译完成后使用 `logout` 退回登录节点。

## MLP 环境搭建

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
pip install triton-2.0.0-1-cp39-cp39-manylinux2014_x86_64.manylinux_2_17_x86_64.whl --no-deps
pip install nvidia_cusparse_cu11-11.7.4.91-py3-none-manylinux1_x86_64.whl --no-deps
pip install nvidia_cuda_cupti_cu11-11.7.101-py3-none-manylinux1_x86_64.whl --no-deps
pip install nvidia_nccl_cu11-2.14.3-py3-none-manylinux1_x86_64.whl --no-deps
pip install nvidia_nvtx_cu11-11.7.91-py3-none-manylinux1_x86_64.whl --no-deps
pip install torch-2.0.1+cu117.with.pypi.cudnn-cp39-cp39-linux_x86_64.whl
pip install mkl_include-2024.2.2-py2.py3-none-manylinux1_x86_64.whl
```

最后，使用 `python3 -c "import torch; print(torch.__version__)"` 验证安装。

### 安装 libtorch

访问 [网址](https://download.pytorch.org/libtorch/cu117/) 下载 `libtorch-shared-with-deps-2.0.1%2Bcu117.zip`。

```
unzip libtorch-shared-with-deps-2.0.1%2Bcu117.zip
```

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

---

## MLP_CPU 环境搭建

首先 `conda create -n mlp_cpu --clone base` 通过克隆 base 环境创建虚拟环境 mlp_cpu。
然后 `conda activate mlp_cpu` 切换虚拟环境。

### 安装 Pytorch

访问 [网址](https://download.pytorch.org/whl/torch/) 下载 Pytorch 安装包 `torch-2.0.1+cpu-cp39-cp39-linux_x86_64.whl`。
之所以不选择 cxx11-abi 版本是因为 `strings /usr/lib64/libstdc++.so.6 | grep GLIBCXX` 发现集群系统的 libstdc++.so.6 版本最高支持到 GLIBCXX_3.4.19，而支持 cxx11-abi 的 PyTorch 通常需要至少 GLIBCXX_3.4.20 或更高版本。

```

```

### 安装 libtorch

访问 [网址](https://download.pytorch.org/libtorch/cpu/) 下载 `libtorch-shared-with-deps-2.0.1%2Bcpu.zip`。

```
unzip libtorch-shared-with-deps-2.0.1%2Bcpu.zip
```

---

## 安装 Lammps

首先下载 Lammps：

```
git clone --depth=1 https://github.com/lammps/lammps
unzip lammps-release.zip
```

然后安装 pair_allegro-main 和 NEP_CPU：

```
git clone --depth 1 https://github.com/mir-group/pair_allegro.git
unzip pair_allegro-main.zip
cd pair_allegro-main/
./patch_lammps.sh ../lammps-release/
```

```
git clone https://github.com/brucefan1983/NEP_CPU.git
unzip NEP_CPU-main.zip
cd NEP_CPU-main/
cp src/* interface/lammps/USER-NEP/
cp interface/lammps/USER-NEP/ ../lammps-release/src/
```

最后，编译安装 LAMMPS：

```
cd ~/soft/lammps-release/
mkdir build
cd build/
cmake ../cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH=$HOME/soft/libtorch -DMKL_INCLUDE_DIR=$HOME/anaconda3/envs/mlp/include -DPKG_PLUGIN=ON -DPKG_NEP=ON
logout
ssh cu02
module load intel/2020.1.217
cd ~/soft/lammps-release/build
make -j$(nproc)

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

