
# CP2K 编译

集群 3 中默认安装了 CP2K 7.1 和 CP2K 9.1。
CP2K 9.1 之后又以每年一到两个更新的速度陆续推出了：

- 2022.1
- 2023.1
- 2023.2
- 2024.1

如果你需要新装一个 CP2K，显然应当安装 2024.1。
特别是，如果你想要在杂化泛函计算时考虑 K 点，那么你的版本号必须 ≥ 2024.1。

CP2K 的安装有多种方法。最简单的方法是直接使用官网预编译版本的二进制可执行文件 `ssmp`。
预编译版本的 CP2K 运行可靠，但由于采用了没有优化过的 `BLAS` 库和 `LAPACK` 库，计算速度较慢。

下面介绍几种不同的自定义安装方法。

## 使用 Apptainer 安装 Docker 打包的镜像

首先在 `Docker Desktop` 中 PULL `apptainer/1.0.0`（因为集群 3 也是这个版本），然后按 `Win + R` 打开 `cmd` 终端，输入：
```sh
docker run -it kaczmarj/apptainer pull docker://cp2k/cp2k:2024.1_openmpi_skylake-avx512_psmp
```
下载 OpenMPI 版本的 `CP2K 2024.1`。将 `SIF` 文件拖入集群，修改提交脚本 `runcp2k.pbs` 即可。

---

## 在虚拟机上使用Docker编译CP2K

首先准备好 Git，执行以下命令：
```sh
git config --global http.postBuffer 524288000
git clone --recursive https://github.com/cp2k/cp2k.git
```
修改 `install_cp2k_toolchain.sh` 中：
```sh
export TARGET_CPU="skylake-avx512"
```
修改 `install_elpa.sh`：
```sh
../configure --prefix="${pkg_install_dir}/${TARGET}/" \
    --libdir="${pkg_install_dir}/${TARGET}/lib" \
    ……
    LIBS="${SCALAPACK_LIBS} $(resolve_string "${MATH_LIBS}" "MPI")" \
    --enable-runtime-threading-support-checks \
    --enable-allow-thread-limiting \
    --without-threading-support-check-during-build \
    > configure.log 2>&1 || tail -n ${LOG_LINES} configure.log
```
修改 `build_dockerhub_images.sh` 中 Dockerfile 的文件名：
```sh
docker build --shm-size=1g --build-arg "GIT_COMMIT_SHA=${SHA}" -f Dockerfile.test_psmp -t "${TAG}" ../../
```
使用 bash 命令运行脚本：
```sh
bash build_dockerhub_images.sh
```
脚本会在虚拟机上生成一个 CP2K 的镜像。使用 `docker images` 获取镜像 ID（例如：913358d296d9），并将其保存为 tar 文件：
```sh
docker save 913358d296d9 -o cp2k-2024.tar
```
将 tar 文件上传至集群，`ssh cu02` 切换到编译节点，加载 apptainer 模块并构建 sandbox：
```sh
module load apptainer/1.0.0
apptainer build --sandbox cp2k-2024 docker-archive://cp2k-2024.tar
```
使用 `apptainer shell -w cp2k-2024` 以写入模式（-w）打开容器，并键入 `nano /.singularity.d/env/91-environment.sh` 编辑环境变量文件如下：
```sh
source /opt/cp2k-toolchain/install/setup
export PATH=$PATH:/opt/cp2k/exe/local
export CP2K_DATA_DIR=/opt/cp2k/data
```
修改完成后，按 `Ctrl + X` 保存并退出。最后，`apptainer build cp2k-2024.sif cp2k-2024` 将容器压缩成 SIF 可执行文件。
这样就完成了 CP2K 的编译和部署，确保在集群上能够充分利用 AVX512 指令集，获得更好的性能。

---

## 从源码编译CP2K的详细步骤

访问 [CP2K官网](https://www.cp2k.org/download) 下载 .tar.bz2 安装包，拖入集群，使用命令 `tar -jxvf cp2k-2023.1.tar.bz2` 解压。
`ssh cu02` 切换到编译节点，然后进入安装目录。

1. 首先打扫干净安装环境：
```sh
make clean
make distclean
rm arch/*
```
2. 然后加载编译环境：
```sh
module load intel/2020.1.217
module load gcc/9.3
```
加载环境后，使用命令 `gcc -v`、`g++ -v`、`gfortran -v` 检查编译器版本。

3. 使用 toolchain 方法安装依赖库：
```sh
cd tools/toolchain/
./install_cp2k_toolchain.sh
```
toolchain 方法在编译时需要联网，`install_cp2k_toolchain.sh` 脚本会自动下载安装需要的依赖库。
`./install_cp2k_toolchain.sh -h` 可以查看帮助文件，其中系统环境变量里已经有的库用 `--with-***=system`；没有安装的库用 `--with-***=install`，脚本会自动联网下载安装；不想安装的库用 `--with-***=no` 即可。
一般来说系统默认的配置方案就是合理的，即什么都不用指定，但是考虑到集群 3 的实际工作情况，在离线编译时，首先 `mkdir build`，然后在 `build/` 中上传已经下载好的对应依赖库的压缩包，然后键入命令（推荐安装 elpa，但是集群 3 安装 elpa 容易失败）：
```sh
nohup ./install_cp2k_toolchain.sh --with-cmake=install --with-openblas=no --with-scalapack=no --with-elpa=no --with-plumed=install --with-sirius=no > install.log 2>&1 &
```
4. 编译 CP2K：
```sh
cp ./install/arch/* ../../arch/
source ./install/setup
cd ../../
make -j 32 ARCH=local VERSION="popt psmp"
```
注意 `make -j` 后面的数字是核数，使用命令 `nproc` 或 `lscpu` 即可查看当前节点 CPU 的核数（登录节点一般是 16 核，编译节点一般是 48 核）。
完成之后使用命令 `make -j 48 ARCH=local VERSION="popt psmp" test` 测试 CP2K。至此，编译和测试完成！

5. 运行 CP2K：
把以下内容加入到 PBS 脚本里：
```sh
module load intel/2020.1.217
module load gcc/9.3
source ~/Softwares/cp2k-2023.1/tools/toolchain/install/setup
export PATH=$PATH:~/Softwares/cp2k-2023.1/exe/local
mpirun -n 24 cp2k.popt cp2k.inp 1>cp2k.out 2>cp2k.err
```
除了使用 PBS 脚本，运行 `cp2k.ssmp -v` 可以查看CP2K的版本、编译时用的库和参数信息；如果要实时查看输出使用命令 `mpirun -np 4 cp2k.popt test.inp | tee test.out` 即可。
