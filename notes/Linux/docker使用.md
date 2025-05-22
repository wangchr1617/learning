## Docker基础操作

```
# 启动docker
systemctl start docker
# 设置开机启动
systemctl enabled docker
# 查看docker信息
docker info

docker run			# 运行容器
docker ps			# 列出正在运行的容器
docker build		# 构建docker镜像
docker exec			# 在容器中执行命令
docker images		# 查看已加载的镜像

# 强制删除镜像
docker rmi -f 293e96aba527

# 启动一个容器
docker run -d ubuntu

# 加载下载好的镜像文件
docker load < docker-oowy-glibc-latest-amd64.tar
```

```
# 运行glibc镜像
docker run -it oowy/glibc /bin/bash
```

## 基于dockerfile配置

首先配置镜像源

```
vim /etc/docker/daemon.json
添加 https://docker.xuanyuan.me

# 重启服务
systemctl daemon-reload
systemctl restart docker
```

联网：/home/zfli-ICME/scripts/auth_client -u by2201136 -p LZFlzf0513

进行配置的dockerfile文件在 **/opt/docker-images** 下。

下面是dockerfile文件的内容

```
FROM ubuntu:24.04

# 避免 systemd 相关问题
ENV DEBIAN_FRONTEND=noninteractive
ENV INITRD=No
ENV container=docker

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

# 安装 Python 3.12
RUN apt install -y python3.12 python3.12-dev python3.12-venv && \
    ln -s /usr/bin/python3.12 /usr/local/bin/python && \
    ln -s /usr/bin/python3.12 /usr/local/bin/python3

RUN python3.12 -m venv /opt/venv
ENV PATH="/opt/venv/bin:$PATH"
RUN /opt/venv/bin/python -m pip install --upgrade pip -i https://pypi.tuna.tsinghua.edu.cn/simple && \
    /opt/venv/bin/pip install ovito -i https://pypi.tuna.tsinghua.edu.cn/simple && \
    /opt/venv/bin/pip install ase -i https://pypi.tuna.tsinghua.edu.cn/simple && \
    /opt/venv/bin/pip install pymatgen -i https://pypi.tuna.tsinghua.edu.cn/simple && \
    /opt/venv/bin/pip install matplotlib -i https://pypi.tuna.tsinghua.edu.cn/simple

# 验证环境
RUN python -c "import ovito; print('OVITO版本:', ovito.version)" && \
    ldd --version | head -n 1

# 设置工作目录
WORKDIR /app

# 启动命令（示例）
CMD ["bash"]
```

在dockefile文件所在目录下运行：`docker build -t ovito_python_3.12 .` 开始构建环境。

构建成功后通过以下命令运行环境：`docker run -it ovito_python_3.12:latest`

更多时候需要挂载自己的目录到容器中的某个目录：`docker run -it -v ~/:/app ovito_python_3.12:latest`

### 通过apptainer在计算节点上配置

```shell
# 在root下先把原来的docker镜像导出为.tar文件，并修改用户权限和所属
docker save -o ovito_python_3.12.tar ovito_python_3.12:latest
cp ovito_python_3.12.tar /home/zfli-ICME/software/docker-images/
chown zfli-ICME:zfli-ICME /home/zfli-ICME/software/docker-images/ovito_python_3.12.tar

# 下面用普通用户权限操作
/home/zfli-ICME/software/docker-images/
chmod 755 ovito_python_3.12.tar
module load apptainer
# 将tar文件转换为沙盒目录
apptainer build --sandbox ovito_sandbox docker-archive:///home/zfli-ICME/software/docker-images/ovito_python_3.12.tar
# 把沙盒目录打包为apptainer镜像文件格式
apptainer build ovito_python.sif ovito_sandbox
# 启动交互式窗口，或者将主机目录挂载到容器内
apptainer shell ovito_python.sif
apptainer run --bind /host/path:/container/path my_image.sif
# 或者执行特定命令（可以在pbs提交脚本中执行）
apptainer exec --bind ~/:/app ~/software/docker-images/ovito_python.sif python /app/scripts/lammps/analysis_trj.py 1600

# 当前目录挂载
apptainer exec --bind $(pwd):/app ovito_python.sif python /app/script.py
```

