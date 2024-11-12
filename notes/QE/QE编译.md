
# QE 编译

## 集群 3 编译 QE 7.0

编译前需要`ssh cu02`切换到编译节点（cu02），并加载编译环境如下：

```sh
module load intel/2020.1.217
module load gcc/9.3
```

执行以下命令查看编译环境是否加载成功：

```sh
icc --version
icpc --version
ifort --version
mpirun --version
echo $MKLROOT
```

输出以下信息即可：

<div align="left">
<img src="./figures/编译_001.png" width = "50%" />
</div>

访问 [GitHub QE 仓库](https://github.com/QEF/q-e/releases) 下载 QE 安装包，以 `qe-7.0-ReleasePack.tgz` 为例。
依次执行以下命令：

```
tar -zvxf qe-7.0-ReleasePack.tgz
cd qe-7.0/
./configure --with-cuda=no --with-scalapack=intel --with-scalapack-qrcp=yes MPIF90=mpiifort FC=ifort FCFLAGS=-O3 CC=icc CFLAGS=-O3 LIBDIRS="${MKLROOT}/lib/intel64/ /opt/intel2020/compilers_and_libraries_2020.1.217/linux/mpi/intel64/lib"
make all -j
```

正常结束的话，`qe-7.0/bin` 目录里有各种 `.x` 结尾的可执行程序，代表安装成功。

如果 QE 在编译时试图从 gitlab.com 下载 devicexlib 库的依赖项，但网络连接失败，无法解析 gitlab.com 的主机地址。
可以手动下载依赖项，并将解压后的文件传输到 QE 源代码目录中的 external/devxlib 文件夹。

编译完成后记得使用命令 `logout` 退出编译节点。

