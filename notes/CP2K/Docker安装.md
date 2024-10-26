
### Docker 安装流程

Docker 是一个开源的应用容器引擎，可以把应用以及依赖包打包到一个轻量级、可移植的容器中，然后发布到 Linux 机器上，实现虚拟化。

容器完全使用沙盒机制，相互之间不会存在任何接口。Docker 可用于开发应用、交付应用、运行应用等场景。

首先打开 Windows `控制面板` `程序和功能` `启用或关闭 Windows 功能`，勾选`打开Hyper-V`、`适用于 Linux 的 Windows 子系统`、`虚拟机平台`三个子功能。

更新成功后重启电脑，以管理员的身份打开 Windows PowerShell，并在终端输入`wsl --update`。
进度条拉满后输入`wsl –install`（这一步可能要翻墙，`-d`参数指定版本，默认是Ubuntu），等待安装完成后填写用户名和密码（密码屏幕不显示），然后再次重启电脑。

访问 Docker 官网下载`Docker Desktop for Windows`安装包，双击执行安装即可。
Docker Desktop 可以通过下载镜像 Images（类似于应用加运行环境的安装包）来 run 一个甚至多个 Containers（类似于应用市场里的应用），数据保存在 Volumes 中。

---

### Docker 的使用
