
# WSL 及 Docker 安装

## WSL 2 安装

WSL（Windows Subsystem for Linux）是微软开发的一项技术，允许用户在 Windows 系统中直接运行完整的 Linux 环境，无需虚拟机。
WSL2 是 WSL 的升级版本，提供了更高的性能和完全的系统调用兼容性。

首先打开 Windows `控制面板 - 程序 - 程序和功能 - 启用或关闭 Windows 功能`，勾选 `打开Hyper-V`、`适用于 Linux 的 Windows 子系统`、`虚拟机平台` 三个子功能。
更新成功后重启电脑确保所有更改生效。

以管理员的身份打开 Windows PowerShell，并在终端输入`wsl --update`。
进度条拉满后输入`wsl –-install`（这一步可能要翻墙，`-d`参数指定版本，默认是Ubuntu），等待安装完成后填写用户名和密码（密码屏幕不显示），然后再次重启电脑。

**再次启动 WSL 时，只需要在终端内键入 `wsl` 即可。**

> 如果 Linux 文件系统不想安装在 C 盘，可以从 [网址](https://cloud-images.ubuntu.com/releases/focal/release/) 下载所需发行版的 `.tar.gz` 文件。
> 假设 `.tar.gz` 文件位于 `D:\Ubuntu` 文件夹下，并且你希望将其导入到 `D:\Ubuntu\Ubuntu-20.04` 路径下，运行：
> 
> ```
> wsl --import MyUbuntu D:\Ubuntu\MyUbuntu D:\Ubuntu\ubuntu-20.04.tar.gz
> ```
> 
> 执行后，WSL 会将 `ubuntu-20.04.tar.gz` 导入到 `D:\Ubuntu\MyUbuntu` 文件夹，并创建该路径下的 Linux 文件系统。
> 可以使用 `wsl --set-default MyUbuntu` 将新的发行版设置为默认。

WSL 安装完成后记得使用 `wsl -l -v` 确保 Ubuntu 发行版已更新为 WSL 2。

首次创建账户成功后可以使用命令 `sudo apt update && sudo apt upgrade -y` 更新软件源并升级所有包；
并使用命令 `sudo apt install -y git curl wget vim` 安装开发常用基础工具。

MobaXterm 提供了对 WSL 的直接支持，无需复杂配置。
打开 MobaXterm，点击顶部工具栏的 "Session" 按钮，在弹出的对话框中选择 "WSL"（在右侧的选项中），点击 OK，MobaXterm 将打开对应的 WSL 终端。

---

## Docker 安装

Docker 是一种流行的容器化平台，它能够简化应用程序的部署和管理。
本文将介绍在 Ubuntu 操作系统上安装 Docker 的步骤，以便我们可以开始使用 Docker 来构建和运行容器化应用程序。

### 安装前的准备工作

如果 Ubuntu 自带 Docker，或者此前安装过低版本的 Docker，可以使用下列命令卸载：
```
sudo apt-get remove docker docker-engine docker.io containerd runc
```

执行以下命令来安装 Docker 依赖：
```
sudo apt install apt-transport-https ca-certificates curl software-properties-common gnupg lsb-release
```

执行以下命令来添加 Docker 阿里镜像的 GPG 密钥：
```
curl -fsSL https://mirrors.aliyun.com/docker-ce/linux/ubuntu/gpg | sudo gpg --dearmor -o /usr/share/keyrings/docker-archive-keyring.gpg
```

添加阿里的 apt 源：
```
echo "deb [arch=$(dpkg --print-architecture) signed-by=/usr/share/keyrings/docker-archive-keyring.gpg] https://mirrors.aliyun.com/docker-ce/linux/ubuntu $(lsb_release -cs) stable" | sudo tee /etc/apt/sources.list.d/docker.list > /dev/null
```

然后更新源：
```
sudo apt update
sudo apt-get update
```

### 安装 Docker

使用下列命令安装最新版本的 Docker：
```
sudo apt install docker-ce docker-ce-cli containerd.io
```

安装完成后，使用 `sudo docker version` 命令查看 Docker 版本，
`sudo systemctl status docker` 命令查看 Docker 运行状态。

然后安装 Docker 命令补全工具：
```
sudo apt-get install bash-completion
sudo curl -L https://raw.githubusercontent.com/docker/docker-ce/master/components/cli/contrib/completion/bash/docker -o /etc/bash_completion.d/docker.sh
source /etc/bash_completion.d/docker.sh
```

当我们安装好了 Docker 之后，需要在 docker 命令前加上 sudo 来执行 docker 命令，例如 `sudo docker ps`。
或者先使用 `sudo -i` 切换至 root 用户，再执行 docker 命令。

如果是在集群上，需要将用户添加到 Docker 用户组，使其无需 sudo 即可运行 Docker 命令。：
```
sudo groupadd docker # 创建 docker 用户组, 注意：Docker 安装后通常会自动创建此组，重复执行可能报错
sudo usermod -aG docker $USER # 将当前用户添加到 docker 组
newgrp docker # 切换当前会话的组到 docker，实际是启动一个新的子 Shell
docker ps -a # 验证用户是否有权限执行 docker 命令，成功列出容器即表示配置正确。
```

最后，我们需要在配置文件中增加国的可用的 Docker hub 镜像。
使用 `sudo vim /etc/docker/daemon.json` 命令编辑配置文件如下：
```
{
  "registry-mirrors": [
    "https://docker.m.daocloud.io", 
	"https://hub-mirror.c.163.com", 
	"https://dockerproxy.com"
  ]
}  
```

然后重启 Docker 服务：
```
sudo service docker restart
```

使用 `docker info` 查看 Docker 信息。

使用 `docker info | grep -A3 "Registry Mirrors"` 查看镜像是否配置成功。

---

## 搜索 Docker 镜像

我们可以从叫做 [Docker hub](https://hub.docker.com/) 的 Docker 官方库获得镜像，或者我们也可以制作自己的镜像。

使用 `docker search` 命令搜索 Docker 镜像，例如： 
```
sudo docker search ubuntu
```

---

## 下载 Docker 镜像

使用 `docker pull` 命令下载 Docker 镜像，例如： 
```
sudo docker pull ubuntu:20.04 # 指定 Ubuntu 镜像版本是 20.04
```

所有已下载的 Docker 镜像都保存在 `/var/lib/docker` 路径下。

使用 `docker images` 命令查看所有已下载的 Docker 镜像，输出如下：
```
wchr@LAPTOP-SPENR6IJ:/mnt/c/Users/wangchangrui$ docker images
REPOSITORY   TAG       IMAGE ID       CREATED       SIZE
ubuntu       20.04     b7bab04fd9aa   6 weeks ago   72.8MB
```

---

## 运行 Docker 容器

基于 Docker 镜像的标签（TAG）或者 ID（Image ID），我们可以使用 `docker run` 命令​​创建并启动指定一个新的容器，例如：
```
sudo docker run -it ubuntu:20.04 /bin/bash # 使用 --name 参数给创建的容器命名
sudo docker run -it b7bab04fd9aa /bin/bash
```
启动容器后，会自动进入容器的 shell（命令行）。

在容器中，可以使用 `exit` 命令终止容器的运行并脱离。
或者按住 `CTRL + P` 然后 `CTRL + Q` 脱离但不关闭容器。

脱离后，使用 `docker ps -a` 命令可以看到当前容器仍在运行（Up）。
```
wchr@LAPTOP-SPENR6IJ:/mnt/c/Users/wangchangrui$ docker ps
CONTAINER ID   IMAGE          COMMAND       CREATED         STATUS         PORTS     NAMES
2105f3ec080a   ubuntu:20.04   "/bin/bash"   3 minutes ago   Up 3 minutes             hardcore_austin
```

当一个新容器被创建后，会自动分配一个唯一的 ID 和 Name。
注意，容器 ID 和 Docker 镜像 ID 是不同的。

我们可以使用 `docker attach` 命令连接到正在运行的容器：
```
sudo docker attach 2105f3ec080a # 基于容器 ID
```
同样按住 `CTRL + P` 然后 `CTRL + Q` 从容器脱离。

使用下面的命令从 Docker 的主机系统中终止正在运行的容器：
```
sudo docker stop 2105f3ec080a # 用空格隔开 ID 可以同时退出多个容器
sudo docker kill 2105f3ec080a # 强制关闭，谨慎使用
```

如果需要启动一个已停止的容器​​（不创建新容器），使用 `docker start` 命令。

---

## 构建自定义 Docker 镜像

我们可以创建自定义的 Docker 镜像，避免不同系统、软件间的环境冲突。

### Dockerfile 构建镜像（以 Ovito 为例）

通过文本文件（Dockerfile）定义镜像的构建步骤，包括基础镜像、依赖安装、文件复制、环境变量配置等。

在一个空目录下，新建一个名为 Dockerfile 文件，并在文件内添加以下内容：
```shell
# 指定基础镜像，用于后续的指令构建。
FROM ubuntu:24.04 

# 在容器内部设置环境变量
ENV DEBIAN_FRONTEND=noninteractive
ENV INITRD=No
ENV container=docker

# 使用 RUN 在镜像中执行 shell 命令。
RUN apt update && apt upgrade -y && \
    apt install -y \
    wget \
    curl \
    git \
    build-essential \
    zlib1g-dev \
    libffi-dev \
    libssl-dev \
    libsqlite3-dev \
    libbz2-dev \
    liblzma-dev \
    libreadline-dev \
    vim

# 安装 Qt6 图形库
RUN apt install -y --no-install-recommends \
    libqt6core6 \
    libqt6gui6 \
    libqt6widgets6 \
    libgl1-mesa-dev

# 安装 Python
RUN apt install -y python3.12 python3.12-dev python3.12-venv && \
    ln -s /usr/bin/python3.12 /usr/local/bin/python && \
    ln -s /usr/bin/python3.12 /usr/local/bin/python3

# 使用轻量级包管理工具 venv
RUN python3.12 -m venv /opt/venv
ENV PATH="/opt/venv/bin:$PATH"
RUN /opt/venv/bin/python -m pip install --upgrade pip -i https://pypi.tuna.tsinghua.edu.cn/simple && \
    /opt/venv/bin/pip install ovito -i https://pypi.tuna.tsinghua.edu.cn/simple && \
    /opt/venv/bin/pip install ase -i https://pypi.tuna.tsinghua.edu.cn/simple && \
    /opt/venv/bin/pip install pymatgen -i https://pypi.tuna.tsinghua.edu.cn/simple && \
    /opt/venv/bin/pip install matplotlib -i https://pypi.tuna.tsinghua.edu.cn/simple

# 验证安装是否成功
RUN python -c "import ovito; print('OVITO版本:', ovito.version)" && ldd --version | head -n 1

# 设置后续指令的工作目录。
WORKDIR /app

# 指定容器创建时的默认命令。
CMD ["bash"]
```

`CMD` 类似于 `RUN` 指令，用于运行程序，但二者运行的时间点不同:
- `CMD` 在docker run 时运行（如果 Dockerfile 中如果存在多个 CMD 指令，仅最后一个生效。）。
- `RUN` 是在 docker build 时运行。

Dockerfile 的指令每执行一次都会在 docker 上新建一层，而过多无意义的层，会造成镜像膨胀过大。
因此，推荐使用 `&&` 符号连接命令，减少创建镜像层数。

在 Dockerfile 文件的存放目录下，执行构建命令：
```
docker build -t ovito_python_3.12:v1 . # （镜像名称:镜像标签）
```

### 直接构建自定义镜像（以 QMCPACK 为例）

首先 `docker pull ghcr.io/qmcpack/ubuntu22-openmpi:latest` 拉取镜像，使用 `docker images` 获取镜像 ID。

启动名为 qmcpack-4.1.0 的临时容器：
```
docker run -it --user root --name qmcpack-4.1.0 ghcr.io/qmcpack/ubuntu22-openmpi:latest /bin/bash
```

克隆 QMCPACK 源代码：
```
git clone https://github.com/QMCPACK/qmcpack.git
```

或者上传下载到本地的 QMCPACK 源代码并解压：
```
docker cp /mnt/e/qmcpack-4.1.0.tar.gz 04cc6d9f0259:/opt/
docker attach 04cc6d9f0259
cd /opt/
tar -xzvf qmcpack-4.1.0.tar.gz
```

编译 QMCPACK
```
cd qmcpack-4.1.0/build
cmake -GNinja -DCMAKE_BUILD_TYPE=RelWithDebInfo -DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpicxx -DQMC_COMPLEX=0 ..  
ninja
ninja install
```

测试 QMCPACK是否编译成功：
```
ctest -VV -R deterministic-unit_test_wavefunction_trialwf
ctest -L deterministic  
```

按住 `CTRL + P` 然后 `CTRL + Q` 脱离容器，
使用命令 `docker commit qmcpack-4.1.0 qmcpack:v4.1.0` 重建镜像。

然后 `docker save -o qmcpack-4.1.0.tar qmcpack:v4.1.0` 保存镜像为 `.tar` 文件。
注意，镜像名称需与重建的镜像一致，生成的 `qmcpack-4.1.0.tar` 文件即为镜像归档文件。

删除镜像命令是 `docker rmi [镜像 ID]`。

---

### 通过 apptainer 在计算节点上配置（以 QMCPACK 为例）

上传 `.tar` 文件到集群，并修改权限：
```
chmod 777 qmcpack-4.1.0.tar 
```

通过 Apptainer 转换为 `.sif` 镜像：
```
module load apptainer
apptainer build --sandbox qmcpack-sandbox docker-archive:///home/changruiwang-ICME/Software/qmcpack/qmcpack-4.1.0.tar
apptainer build qmcpack.sif qmcpack-sandbox
```

由于已通过 Docker 编译 QMCPACK 并固化到镜像中​​，因此​​无需挂载。
​直接运行 `apptainer exec qmcpack.sif qmcpack <参数>` 即可。
