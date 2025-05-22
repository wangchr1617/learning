
# Linux 中的 Module 工具使用指南

Module 是一个专门管理环境变量的工具，通常用于软件或运行库等设备有多个版本时，需要分别配置这些环境变量。
它依赖于 TCL 工具，可以通过二进制编译安装或使用 yum/apt 快速安装。
---

## 1. Module 的安装和初始化

首先，需要安装 TCL 工具：
```
wget https://cfhcable.dl.sourceforge.net/project/tcl/Tcl/8.5.9/tcl8.5.9-src.tar.gz
tar -zxvf tcl8.5.9-src.tar.gz
cd tcl8.5.9/unix
./configure --prefix=/usr/local/tools/tcl
make
make install
```
然后，安装 module 工具：
```
wget https://newcontinuum.dl.sourceforge.net/project/modules/Modules/modules-4.2.4/modules-4.2.4.tar.gz
tar -zxvf modules-4.2.4.tar.gz
cd modules-4.2.4
./configure --prefix=/usr/local/tools/modules --with-tcl-lib=/usr/local/tools/tcl/lib --with-tcl-inc=/usr/local/tools/tcl/include
make
make install
```
安装完成后，需要初始化 module 工具，可以通过 source 命令加载相应的初始化文件：
```
source /usr/local/tools/modules/init/profile.sh
```
或者将其添加到 `/etc/profile` 文件中，以便系统重启后自动加载：
```
ln -s /usr/local/tools/modules/init/profile.sh /etc/profile.d/module.sh
```
---

## Module 的使用

module 工具依托于 MODULEPATH 环境变量来查找配置信息目录。设置好目录结构和环境变量后，只需设置 MODULEPATH，module 工具就会自动查找路径下的所有配置信息。例如：
```
export MODULEPATH=/opt/modulefiles
```
常用的 module 命令包括：
```
module avail # 显示可用的模块
module list # 显示已加载的模块
module load gcc/4.8.4 # 加载模块
module unload gcc/4.8.4 # 卸载模块
module purge # 取消所有加载的工具
module show gcc/4.8.4 # 查看相应配置信息
```

---

## Module 配置文件编写

以 `/opt/modulefiles/gcc/9.3` 为例，
```
#%Module1.0
conflict     GCC-9.3
module-whatis " modulefile about gcc@9.3.0 "
set topdir "/opt/gcc/gcc-9.3.0"
prepend-path MANPATH "/opt/gcc/gcc-9.3.0/share/man"
prepend-path PATH "${topdir}/bin"
prepend-path LIBRARY_PATH "${topdir}/lib"
prepend-path LIBRARY_PATH "${topdir}/lib64"
prepend-path LD_LIBRARY_PATH "${topdir}/lib"
prepend-path LD_LIBRARY_PATH "${topdir}/lib64"
prepend-path CPATH "${topdir}/include"
```

Module 工具中，`setenv`、`append-path` 和 `prepend-path` 是用于操作环境变量的关键命令，它们在 Modulefile 中有不同的用途和效果：
```
# 设置软件根目录（非路径变量）
setenv MYAPP_HOME "/opt/myapp"

# 添加可执行目录到 PATH 开头（优先使用）
prepend-path PATH "$env(MYAPP_HOME)/bin"

# 添加库目录到 LD_LIBRARY_PATH 末尾（备选）
append-path LD_LIBRARY_PATH "$env(MYAPP_HOME)/lib"
```

如果我想指定 gcc/9.3 为默认版本，可以在 `modulefile` 同级目录下创建 `.version` 文件如下所示：
```
#%Module1.0
set ModulesVersion "9.3"
```
