
# WSL 及 Docker 安装

## WSL 2 安装

首先打开 Windows `控制面板 - 程序 - 程序和功能 - 启用或关闭 Windows 功能`，勾选`打开Hyper-V`、`适用于 Linux 的 Windows 子系统`、`虚拟机平台`三个子功能。

更新成功后重启电脑，以管理员的身份打开 Windows PowerShell，并在终端输入`wsl --update`。
进度条拉满后输入`wsl –install`（这一步可能要翻墙，`-d`参数指定版本，默认是Ubuntu），等待安装完成后填写用户名和密码（密码屏幕不显示），然后再次重启电脑。

> 如果 Linux 文件系统不想安装在 C 盘，可以从 [网址](https://cloud-images.ubuntu.com/releases/focal/release/) 下载所需发行版的 `.tar.gz` 文件。
> 假设 `.tar.gz` 文件位于 `D:\Ubuntu` 文件夹下，并且你希望将其导入到 `D:\Ubuntu\Ubuntu-20.04` 路径下，运行：
> 
> ```
> wsl --import MyUbuntu D:\Ubuntu\MyUbuntu D:\Ubuntu\ubuntu-20.04.tar.gz
> ```
> 
> 执行后，WSL 会将 `ubuntu-20.04.tar.gz` 导入到 `D:\Ubuntu\MyUbuntu` 文件夹，并创建该路径下的 Linux 文件系统。
> 可以使用 `wsl --set-default MyUbuntu` 将新的发行版设置为默认。

WSL 安装完成后记得使用 `wsl -l -v` 确保 Ubuntu 发行版已更新为 WSL 2。然后运行 `wsl` 启动 WSL 2。

MobaXterm 提供了对 WSL 的直接支持，无需复杂配置。
打开 MobaXterm，点击顶部工具栏的 "Session" 按钮，在弹出的对话框中选择 "WSL"（在右侧的选项中），点击 OK，MobaXterm 将打开对应的 WSL 终端。
---

## Docker 安装

Docker 是一个开源的应用容器引擎，可以把应用以及依赖包打包到一个轻量级、可移植的容器中，然后发布到 Linux 机器上，实现虚拟化。
容器完全使用沙盒机制，相互之间不会存在任何接口。Docker 可用于开发应用、交付应用、运行应用等场景。

访问 Docker 官网下载 `Docker Desktop for Windows` 安装包，双击执行安装即可。
Docker Desktop 可以通过下载镜像 Images（类似于应用加运行环境的安装包）来 run 一个甚至多个 Containers（类似于应用市场里的应用），**数据保存在 Volumes 中**。
安装并启动 Docker Desktop 后，运行 `docker --version` 检查 Docker 版本。