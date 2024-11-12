
# GROMACS 编译

首先编译快速傅立叶变换库 fftw。
访问 [下载链接](http://www.fftw.org/fftw-3.3.8.tar.gz) 下载压缩包，拖入集群解压，进入此目录后运行：

```
module load gcc/9.3
./configure CC=gcc --prefix=/home/changruiwang-ICME/Software/fftw338 --enable-shared --enable-float --enable-sse --enable-sse2 --enable-avx --enable-avx2 --enable-avx512
```

然后 `make -j install` 开始编译。

访问 [GROMACS 官网](http://ftp.gromacs.org/pub/gromacs/gromacs-2018.8.tar.gz) 下载 GROMACS 2018.8 压缩包，拖入集群解压，进入此目录后运行：

```
module load intel/2020.1.217
mkdir build
cd build
export CMAKE_PREFIX_PATH=/home/changruiwang-ICME/Software/fftw338
chmod u+w /home/changruiwang-ICME/Software/gromacs
cmake .. -DCMAKE_INSTALL_PREFIX=/home/changruiwang-ICME/Software/gromacs/gmx2018.8 -DGMX_MPI=ON -DGMX_SIMD=AVX2_256
make -j
make install
```

即可完成 GROMACS 编译。
