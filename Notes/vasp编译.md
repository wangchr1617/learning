# VASP 编译

注意，不同版本的 VASP 编译过程不同，不同集群的编译环境也有所不同。

---

### 集群3 编译 VASP 6.3.0

编译前需要`ssh cu02`切换到编译节点（cu02），并加载编译环境如下：
```sh
module load intel/2020.1.217
module load gcc/9.3
```
下载 VASP 安装包并解压。`cd vasp.6.3.0/`进入安装目录，复制 `makefile.include` 文件：
```sh
cp arch/makefile.include.intel ./makefile.include
```
如果想要使用 Intel 的 MKL 数学库，需要修改 `makefile.include` 文件中的一行：
```sh
FCL+=-qmkl=sequential  # 原行
FCL+=-mkl=sequential  # 修改后
```

---

### 集群2 编译 VASP 5.4.4

编译前需要`ssh cu10`切换到编译节点（cu10），并加载编译环境如下：
```
source /opt/intel/composer_xe_2015/bin/compilervars.sh intel64
source /opt/intel/mkl/bin/intel64/mklvars_intel64.sh
source /opt/intel/impi/5.0.2.044/bin64/mpivars.sh
```
集群2编译 VASP 5.4.4 时用默认的 gcc 即可。

下载 VASP 安装包并解压。`cd vasp.5.4.4/`进入安装目录，复制 `makefile.include` 文件：
```sh
cp arch/makefile.include.linux_intel ./makefile.include
```

---

### 集群2 注意事项

1. 编译Intel fftw3库

**[2024-07-20] 集群2已经不再需要手动编译Intel fftw3库了**

集群2使用 Intel 2015 时，由于没有编译 fftw3xf 库，需要拷贝 `fftw3xf` 文件夹到自定义文件夹：
```sh
cp /opt/intel/mkl/interfaces/fftw3xf .
```
然后在该文件夹下 `make libintel64` 编译 `libfftw3xf_intel.a` 文件，并将生成的 `.a` 文件路径添加到 `makefile.include` 文件的 `OBJECTS` 行。

2. `libmkl_intel_lp64.so: cannot open shared object file`

通过在提交脚本中设置环境变量 `LD_LIBRARY_PATH` 来解决：
```sh
export LD_LIBRARY_PATH=/opt/intel/mkl/lib/intel64:$LD_LIBRARY_PATH
```

3. `/lib64/libc.so.6: version GLIBC_2.14' not found`
同上，设置环境变量：
```
export LD_LIBRARY_PATH=/opt/glibc-2.14/lib:$LD_LIBRARY_PATH:
```
或在上一步的环境变量后追加一个环境变量：
```
export LD_LIBRARY_PATH=/opt/intel/mkl/lib/intel64:$LD_LIBRARY_PATH:/opt/glibc-2.14/lib/
```
通过 `echo $LD_LIBRARY_PATH` 检查环境变量设置是否正确。


---

### 固定基矢

结构优化过程中固定晶格基矢有两种方法：

**方法一**

这是网上流传最广的补丁方法，使用输入文件 `OPTCELL` 固定。根据 [这篇博文](https://blog.shishiruqi.com//2019/05/05/constr/) 的内容重写 `./constr_cell_relax.F` 文件。`OPTCELL` 的范例如下：
```
100
110
000
```

**方法二**

打开 [VASP_OPT_AXIS 的 GitHub 页面](https://github.com/Chengcheng-Xiao/VASP_OPT_AXIS)，下载 `.zip` 补丁文件。
解压后 `cd VASP_OPT_AXIS-master/fixing_stress_tensor/` 。
在这里可以看到两个补丁文件，其中 `stress_relax.patch` 等同于前述方法一；而 `stress_relax_finner.patch` 通过在 `INCAR` 中设置 `IOPTCELL` 参数来固定晶格基矢。
将 `stress_relax_finner.patch` 文件复制到 `vasp.6.3.0/` 目录，输入以下命令打补丁：
```sh
patch -p0 < stress_relax_finner.patch
```
输出 `succeeded` 表示补丁已经打好。

---

### 集群3 VASP 6.3.0 添加 VTST 功能
为实现过渡态结构的搜索，需要把 VTST 的一些功能添加到 VASP 中。
访问 [VTST 官网](https://theory.cm.utexas.edu/vtsttools/download.html) 下载安装包（例如`vtstcode-198.tgz`）并拖入集群解压。

```sh
cp -r vasp.6.3.0/ vasp-vtst/
cp vasp-vtst/src/chain.F vasp-vtst/src/chain.F.copy
cp vtstcode6.3/* vasp-vtst/src/
cd vasp-vtst/
```
编辑 `src/main.F` 文件，替换以下内容
```fortran
CALL CHAIN_FORCE(T_INFO%NIONS,DYN%POSION,TOTEN,TIFOR, &
     LATT_CUR%A,LATT_CUR%B,IO%IU6)
```
为
```fortran
CALL CHAIN_FORCE(T_INFO%NIONS,DYN%POSION,TOTEN,TIFOR, &
     TSIF,LATT_CUR%A,LATT_CUR%B,IO%IU6)
```
然后替换以下内容
```fortran
IF (LCHAIN) CALL chain_init( T_INFO, IO)
```
为
```
CALL chain_init( T_INFO, IO)
```
编辑 `src/.objects` 文件，在 `chain.o` 之前插入以下内容：
```makefile
bfgs.o dynmat.o instanton.o lbfgs.o sd.o cg.o dimer.o bbm.o \
fire.o lanczos.o neb.o qm.o \
pyamff_fortran/*.o ml_pyamff.o opt.o \
```
编辑 `src/makefile` 文件，找到 `LIB` 变量，修改如下：
```makefile
LIB= lib parser pyamff_fortran
```
然后找到 `dependencies`，修改如下：
```makefile
dependencies: sources libs
```

---

### 集群2 VASP 5.4.4 添加 VTST 功能
为实现过渡态结构的搜索，需要把 VTST 的一些功能添加到 VASP 中。
访问 [VTST 官网](https://theory.cm.utexas.edu/vtsttools/download.html) 下载安装包（例如`vtstcode-198.tgz`）并拖入集群解压。

```sh
cp -r vasp.5.4.4/ vasp-vtst/
cp vasp-vtst/src/chain.F vasp-vtst/src/chain.F.copy
cp vtstcode5/* vasp-vtst/src/
cd vasp-vtst/
```
编辑 `src/main.F` 文件，替换以下内容
```fortran
CALL CHAIN_FORCE(T_INFO%NIONS,DYN%POSION,TOTEN,TIFOR, &
     LATT_CUR%A,LATT_CUR%B,IO%IU6)
```
为
```fortran
CALL CHAIN_FORCE(T_INFO%NIONS,DYN%POSION,TOTEN,TIFOR, &
     TSIF,LATT_CUR%A,LATT_CUR%B,IO%IU6)
```
编辑 `src/.objects` 文件，在 `chain.o` 之前插入以下内容：
```makefile
bfgs.o dynmat.o instanton.o lbfgs.o sd.o cg.o dimer.o bbm.o \
fire.o lanczos.o neb.o qm.o opt.o \ 
```

---

### 添加 VASPsol 功能
访问[网址](https://github.com/henniggroup/VASPsol)下载安装包并拖入集群解压。
```sh
cp solvation.F ../vasp-vtst/src/solvation.F
```
并在`makefile.include`中的`CPP_OPTIONS=`后面添加`-Dsol_compat \`即可。

---

### 编译

然后使用以下命令开始编译：
```sh
make -j all
# 或
make DEPS=1 -j36 <target> # 这里的 all 和 <target> 包括 std、gam、ncl
```

在`makefile.include`中指定`OFLAG = -O3`可以提高编译性能。

编译完成后记得`logout`退出编译节点：

---