
# GPUMD 编译指南

---

## 1. 下载和编译 GPUMD

1. 访问 [GPUMD GitHub 页面](https://github.com/brucefan1983/) 并下载 `.tar.gz` 安装包。
2. 将安装包上传到集群并解压：
    ```bash
    tar -xzvf <package_name>.tar.gz
    ```
3. 进入 `src` 目录并编译：
    ```bash
    cd <package_name>/src/
    make
    ```
4. 编译完成后，将产生两个可执行文件：`gpumd` 和 `nep`。

5. 环境变量设置

如果编译失败，大概率是因为 `~/.bashrc` 里 `Intel` 和 `GCC` 相关的环境变量设置有问题。
我的参考设置如下：
```
#### ---------------- Intel ------------------- ####
export INTEL_MPI_HOME=/opt/intel2020/compilers_and_libraries_2020.1.217/linux/mpi
export PATH=$INTEL_MPI_HOME/intel64/bin:$PATH
export MANPATH=$INTEL_MPI_HOME/man:$MANPATH
export LD_LIBRARY_PATH=$INTEL_MPI_HOME/intel64/lib:$LD_LIBRARY_PATH
# source /opt/intel2020/compilers_and_libraries_2020.1.217/linux/bin/compilervars.sh intel64
# source /opt/intel2020/mkl/bin/mklvars.sh intel64
# source /opt/intel2020/intelpython3/bin/mpivars.sh
###----------------------------------------------

#### ----------------- GCC -------------------- ####
export PATH=/opt/gcc/gcc-9.3.0/bin:${PATH}
export LIBRARY_PATH=/opt/gcc/gcc-9.3.0/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/opt/gcc/gcc-9.3.0/lib64:$LD_LIBRARY_PATH
export C_INCLUDE_PATH=/opt/gcc/gcc-9.3.0/include:$C_INCLUDE_PATH
export CPLUS_INCLUDE_PATH=/opt/gcc/gcc-9.3.0/inclde:$CPLUS_INCLUDE_PATH
###----------------------------------------------
```

---

## 2. GPUMD 提交脚本

以下是一个示例的 GPUMD 提交脚本：

```bash
#PBS -N gpumd
#PBS -l nodes=1:ppn=32
#PBS -l walltime=144:00:00
#PBS -q gpu
#PBS -S /bin/bash
#PBS -V

NP=`cat $PBS_NODEFILE | wc -l`
NN=`cat $PBS_NODEFILE | sort | uniq | tee /tmp/nodes.$$ | wc -l`
cat $PBS_NODEFILE > /tmp/nodefile.$$

cd $PBS_O_WORKDIR
ulimit -s unlimited
export CUDA_VISIBLE_DEVICES=0,1 # 指定 CUDA 程序可以访问的 GPU 设备

EXEC_gpumd=/home/changruiwang-ICME/Software/GPUMD-3.9.1/src/gpumd
EXEC_nep=/home/changruiwang-ICME/Software/GPUMD-3.9.1/src/nep

$EXEC_gpumd > output # 运行 GPUMD，训练时换成 $EXEC_nep 即可
```

---

## 3. 编译带 NEP 插件的 LAMMPS

### 3.1 加载编译器

首先，加载需要的编译器模块：

```bash
ssh cu02
module load intel/2020.1.217
module load gcc/9.3
```

### 3.2 下载并解压安装包

访问 [NEP_CPU GitHub 页面](https://github.com/brucefan1983/NEP_CPU) 下载 `.zip` 安装包。
访问 [lammps GitHub 页面](https://github.com/lammps/lammps) 下载 `.tar.gz` 安装包，最好是 stable 版本。

```bash
unzip NEP_CPU-main.zip
tar -xzvf lammps-stable.tar.gz
```

### 3.3 准备 NEP 插件

```bash
cp NEP_CPU-main/src/* NEP_CPU-main/interface/lammps/USER-NEP/
cp NEP_CPU-main/interface/lammps/USER-NEP/ lammps-2Aug2023/src/
cd lammps-2Aug2023/src
make yes-USER-NEP
```

### 3.4 编译 LAMMPS

```bash
mkdir ../build
cd ../build
cmake -D PKG_BODY=ON -D PKG_FEP=ON -D PKG_EXTRA-COMPYTE=ON -D PKG_KSPACE=ON -D PKG_MISC=ON -D PKG_MOLECULE=ON -D PKG_OPT=ON -D PKG_REAXFF=ON -D PKG_REPLICA=ON -D PKG_RIGID=ON -D BUILD_SHARED_LIBS=yes ../cmake
cmake --build . -j 4
```

编译完成后，使用 `logout` 命令退回到登录节点。

---

## 4. LAMMPS 提交脚本

以下是一个示例的 LAMMPS 提交脚本：

```bash
#PBS -N lammps
#PBS -l nodes=3:ppn=24
#PBS -l walltime=600:00:00
#PBS -q batch
#PBS -V
#PBS -S /bin/bash

module load intel/2020.1.217

NP=`cat $PBS_NODEFILE | wc -l`
NN=`cat $PBS_NODEFILE | sort | uniq | tee /tmp/nodes.$$ | wc -l`
cat $PBS_NODEFILE > /tmp/nodefile.$$

cd $PBS_O_WORKDIR
ulimit -s unlimited
export OMP_NUM_THREADS=3

EXEC=/home/changruiwang-ICME/Software/lammps_nep/build/lmp
mpirun -machinefile $PBS_NODEFILE -np $NP $EXEC -in in* > output
```
