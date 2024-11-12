
# OPENMX 编译

下载资源：
```
wget https://www.openmx-square.org/openmx3.9.tar.gz
tar -xvf openmx3.9.tar.gz 
cd openmx3.9/source
wget https://www.openmx-square.org/bugfixed/21Oct17/patch3.9.9.tar.gz
tar -xvf patch3.9.9.tar.gz
mv kpoint.in ../work/
```

修改 `openmx3.9/source/makefile` 第 8 到 11 行如下：
```
MKLROOT = /opt/intel2020/compilers_and_libraries_2020.1.217/linux/mkl
CC = mpiicc -O3 -xHOST -ip -no-prec-div -qopenmp -I${MKLROOT}/include -I${MKLROOT}/include/fftw
FC = mpiifort -O3 -xHOST -ip -no-prec-div -qopenmp
LIB= -L${MKLROOT}/lib/intel64 -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lifcore -lmkl_blacs_intelmpi_lp64 -liomp5 -lpthread -lm -ldl 
```

开始编译：

```
module load intel/2020.1.217
make clean
make all -j 4 
make install
```
