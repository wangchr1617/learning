# Git 入门教程

## Git 的安装
Git 是 Linus Torvalds 为了帮助管理 Linux 内核开发而开发的一个开放源码的分布式版本控制系统。
访问 [Git 官网](https://git-scm.com/) 下载 Git 安装包。

以 Windows 系统为例，以管理员身份运行 .exe 客户端进入安装界面。
除了修改安装路径，其余的选项保持默认即可。

安装完成后，打开 git bash 输入 `git --version` 验证是否安装成功。

想了解 Git 的各式工具该怎么用，可以阅读它们的使用帮助，方法有三：
```
git help <verb>
git <verb> --help
man git-<verb>
```
比如，要学习 config 命令可以怎么用，运行：
```
git help config
```

## 注册账号
通过访问 [Github 首页](https://github.com/) 或者 [Gitee 首页](https://gitee.com/) 注册个人账号。

## 配置 ssh 公钥
目前 Github 和 Gitee 支持使用 HTTPS 协议和 ssh 协议进行代码的推送/拉取。

使用 https 协议对初学者来说会比较方便，复制 https url 然后到 git Bash 里面直接用 clone 命令克隆到本地就好了，但是每次 fetch 和 push 代码都需要输入账号和密码，这也是 https 协议的麻烦之处。

Github 和 Gitee 均提供了基于 SSH 协议的 Git 服务，在使用 SSH 协议访问仓库之前，需要先配置好账户/仓库的 SSH 公钥。
首先，使用 `git config --global --list` 查看全局设置。
如果已经配置 Git ，需要先清除全局设置：
```
git config --global --unset user.name
git config --global --unset user.email
```

然后 `cd ~/.ssh` 进入 .ssh 文件夹下使用命令 `ssh-keygen -t ed25519 -C "xxxxx@xxxxx.com"` 生成 ssh key 。
注意：这里的 xxxxx@xxxxx.com 只是生成的 ssh key 的名称，并不约束或要求具体命名为某个邮箱。但为了便于辨识，建议使用注册账号的邮箱命名 ssh key 。
按照提示完成三次回车，即可生成 ssh key。

如果需要同时配置 Github 和 Gitee 的 Git ssh key ，不妨在第一个提示后设置 Github 的 ssh key 为 id_github；Gitee 的 ssh key 为 id_gitee 。
完成后会在 ~/.ssh 目录下生成以下四个文件：
- id_github
- id_github.pub
- id_gitee
- id_gitee.pub

复制 id_github.pub 和 id_gitee.pub 的文件内容，分别在 Github 和 Gitee 网页端添加 SSH 公钥。
以 Gitee 为例，通过仓库主页 「管理」->「部署公钥管理」->「添加部署公钥」 ，添加生成的 public key 到仓库中。

在 ~/.ssh 文件夹中创建 config 文件并添加以下内容以区分两个 ssh key ：
```
# Github 配置
Host github.com
    User git
    IdentityFile ~/.ssh/id_github
    HostName github.com
    PreferredAuthentications publickey
    
# Gitee 配置
Host gitee.com
    User git
    IdentityFile ~/.ssh/id_gitee
    HostName gitee.com
    PreferredAuthentications publickey
```

`vi ~/.gitconfig` 编辑 Git 全局配置文件，然后在其中添加以下内容，根据主机名设置不同的 user.name 和 user.email：
```
[user]
    name = Default Username
    email = default@example.com

# 针对 github.com 的特殊配置
[includeIf "gitdir:~/path/to/github-repo/"]
    path = ~/.gitconfig-github

# 针对 gitee.com 的特殊配置
[includeIf "gitdir:~/path/to/gitee-repo/"]
    path = ~/.gitconfig-gitee
```
注意，includeIf 中的 gitdir 路径可以使用通配符。
例如，如果你的 GitHub 仓库都存放在 ~/projects/github/ 目录下，可以写成：
```
[includeIf "gitdir:~/projects/github/"]
    path = ~/.gitconfig-github
```

最后，为 GitHub 和 Gitee 创建两个单独的 Git 配置文件 ~/.gitconfig-github 和 ~/.gitconfig-gitee：
```
# ~/.gitconfig-github
[user]
    name = GitHub Username
    email = github@example.com
```

```
# ~/.gitconfig-gitee
[user]
    name = Gitee Username
    email = gitee@example.com
```

创建 SSH Key 并添加到 GitHub 和 Gitee 后，可以在终端使用以下命令来测试连接：
```
ssh -T git@github.com
ssh -T git@gitee.com
```

首次使用需要确认并添加主机到本机 SSH 可信列表。若返回 `Hi XXX! You've successfully authenticated, but Gitee.com does not provide shell access.` 内容，则证明添加成功。
添加成功后，就可以使用 SSH 协议对仓库进行操作了。

注意，如果是从Windows系统提交，不妨使用
```
git config --global core.autocrlf input
```
避免由于回车（CR）和换行（LF）导致的提交错误。

## 新建仓库
在注册完成并成功登录 Gitee 账号后，用户可以开始创建自己的第一个仓库。

在新建仓库页面填写仓库信息。仓库相关概念说明如下：

- **仓库名称**： 仓库的名称，用于仓库命名。
- **归属**：仓库归属账户，可以是个人账号/组织/企业中的一种，创建成功后该账户默认为仓库的拥有者（管理员）。
- **路径**：仓库的 Git 访问路径，由用户个性地址 + 仓库路径名称组成。创建仓库后用户将通过该路径访问仓库。
- **仓库介绍**：仓库的简单介绍。
- **是否开源**：设置仓库是否为公开仓库，公开仓库对所有人可见，私有仓库仅限仓库成员可见。
- **选择语言**：仓库主要开发用的编程语言。
- **添加.gitignore**：系统默认提供的 Git 忽略提交的文件模板，设置 .gitignore 后将默认忽略指定目录/文件到仓库。
- **添加开源许可证**：如果仓库为公开仓库，可以添加设置仓库的开源协议，作为对当前项目仓库和衍生项目仓库许可约束，开源许可证决定了该开源项目是否对商业友好。
- **Readme**：项目仓库自述文档，通常包含有软件的描述或使用的注意事项。
- **使用模板文件初始化仓库**：使用 Issue 或 Pull Request 文件模板初始化仓库。

点击「创建」，即可在 Gitee 上创建你的仓库。

## 提交代码
**方法1、先将仓库 clone 到本地，修改后再 push 到 Gitee 的仓库**
```
git clone XXX
```
如果需要 clone 指定的分支，可以在仓库地址后面加上分支名，例如：
```
git clone XXX dev
```

修改代码后，可以使用 Git Bash 执行下面命令：
```
git add . # 提交当前仓库的所有改动到 Git 暂存区
git status -s # 查看仓库当前文件提交状态（A：提交成功；AM：文件在添加到缓存之后又有改动）
git commit -s -m "1.0.0" # 从 Git 的暂存区提交版本到仓库，参数 -m 后为当次提交的备注信息
git push origin master # 将本地的 Git 仓库信息推送上传到服务器，参数 origin 后面是提交到的具体分支
git log # 查看 Git 提交的日志
```
撰写提交说明时，如果一行不够，可以只执行 `git commit` ，就会跳出文本编译器，编写多行提交说明。

**方法2、本地初始化一个仓库，设置远程仓库地址后再做push**

和方法1的差别，在于先在本地创建仓库。
```
git init # 初始化新仓库，在当前目录下会出现一个名为 .git 的目录
git remote add origin https://gitee.com/用户个性地址/HelloGitee.git
git remote -v # 查看当前仓库对应的远程仓库地址
```
通过将一个远程仓库添加到本地的仓库中，完成了版本的一次初始化。

接下去，进入你已经初始化好的或者克隆仓库的目录,然后执行：
```
git pull origin master # origin 是仓库名；master 是分支名
```
修改/添加文件，否则与原文件相比就没有变动。
```
git add .
git commit -s -m "第一次提交"
git push origin master
```
在新建仓库时，如果在 Gitee 平台仓库上已经存在 readme 或其他文件，在提交时可能会存在冲突，这时用户需要选择的是保留线上的文件或者舍弃线上的文件，如果您舍弃线上的文件，则在推送时选择强制推送。

强制推送示例：
```
git push origin master -f
```
如果您选择保留线上的 readme 文件,则需要先执行：
```
git pull origin master
```

## 分支管理
`git branch` 是 Git 中管理分支的核心命令之一，它允许你查看、创建、重命名和删除分支。
分支是 Git 中用于管理不同开发线路的工具，便于进行功能开发、错误修复、多人协作等工作。

使用`git branch` 查看分支列表，当前所在的分支会以 * 号标记。
使用 `git branch <branch-name>` 创建新分支，并使用 `git checkout <branch-name>` 切换分支。
通过选项控制如 `git branch -option <new-branch-name>` 可以实现：
-m 重命名当前所在分支；-d 删除本地分支（前提是该分支已合并到当前分支）。

分支合并需要先切换到目标分支（如 main），然后将指定分支（如 dev）合并到当前分支：
```
git checkout main
git merge dev
```

## 撤销命令
撤销命令使用是非常频繁的，因为某些原因，我们不再需要我们的改动或者新的改动有点问题，我们需要回退到某个版本，
这时就需要用到撤销命令，或者说这个应该翻译成重置更加恰当。具体命令如下:
```
git reset --hard 版本号 # 版本号使用 git log 查看
```

## 更新命令
`git fetch` 是 Git 中用于从远程仓库获取更新的命令。
它不会对当前的工作目录做任何更改，而是将远程仓库的更改（如新提交、新分支、变更的引用等）下载到本地。
这使得 `git fetch` 成为一个安全的获取更新的方式，因为它不会直接将远程的更改合并到当前分支，而只是更新本地对远程仓库的引用（如 origin/main）。

`git fetch` 与 `git pull` 的区别在于：
`git pull` 会直接将远程的更改合并到当前分支，
而 `git fetch` 只会下载更新而不进行合并操作。
你可以在 fetch 之后查看、比对差异，然后决定如何处理这些更新。

## 使用 Git Tag 标记版本

以 `main` 分支为例，要将当前 `main` 分支的状态发布为 `1.0.0` 版本，通常使用 Git Tag 来标记版本号。
Tag 是 Git 中用于标识某个提交的特定点，可以用来表示版本。
相较于分支，Tag 更适合用来标记软件发布版本。

首先，检查当前的分支和状态，确保当前在 `main` 分支，并且所有的修改都已经提交。
```
git checkout main
git status
```

使用 `git tag` 命令为当前 `main` 分支的状态创建 `1.0.0` 版本 Tag。
```
git tag -a v1.0.0 -m "Release version 1.0.0"
```
这会创建一个名为 `v1.0.0` 的 Tag，并附上版本说明 "Release version 1.0.0"。

将刚才创建的 Tag 推送到远程仓库，以便在 GitHub 上查看。
```
git push origin v1.0.0
```

使用以下命令查看本地和远程的 Tag 列表，确认 Tag 已成功创建并推送。
```
git tag        # 查看本地 Tag 列表
git ls-remote --tags origin  # 查看远程 Tag 列表
```
