# VSCode + WSL + Claude Code 开发环境配置

本文介绍如何在 Windows 系统中搭建一套完整的 AI 辅助开发环境：从 WSL 和 VSCode 的安装，到 Claude Code 的配置与使用。

---

## 环境信息

| 项目 | 值 |
|---|---|
| 操作系统 | Windows 11 + WSL 2 (Ubuntu) |
| Node.js | v22.22.0（通过 nvm 管理） |
| Claude Code | VSCode 扩展 + npm 全局安装 |
| ccswitch | CC-Switch-v3.16.2（Windows MSI 安装） |
| API 后端 | DeepSeek（`api.deepseek.com`） |
| 代理 | Clash（端口 7890） |

---

## 一、Windows 中安装 VSCode 和 WSL

### 1. 安装 VSCode

访问 [VSCode 官网](https://code.visualstudio.com/) 下载安装包，安装时勾选 **"添加到 PATH"**，以便在终端中直接使用 `code` 命令。

### 2. 安装 WSL 2

WSL 2 的安装详见 [WSL 及 Docker 安装](https://github.com/wangchr1617/learning/blob/main/notes/%E6%96%B0%E7%94%9F%E5%85%A5%E9%97%A8/WSL%E5%8F%8ADocker%E5%AE%89%E8%A3%85.md)，此处仅列关键步骤：

1. 控制面板 → 启用或关闭 Windows 功能 → 勾选 **Hyper-V**、**适用于 Linux 的 Windows 子系统**、**虚拟机平台**
2. 管理员 PowerShell：`wsl --update` → `wsl --install`
3. 重启电脑，WSL 终端中设置用户名和密码

安装完成后用 `wsl -l -v` 确认版本为 WSL 2：

```
C:\Users\xxx> wsl -l -v
  NAME      STATE           VERSION
* Ubuntu    Running         2
```

---

## 二、Windows 中 Clash 的代理设置

WSL 中的很多操作（如 npm 安装、git clone）需要科学上网。本节说明 Windows 端 Clash 的配置要点。

### 1. 开启 Clash 并允许局域网

- 启动 Clash，勾选 **"允许局域网"（Allow LAN）**
- 查看端口号：默认 **7890**（HTTP）和 **7891**（SOCKS5）

### 2. 获取宿主机 IP

在 Windows 的 cmd 或 PowerShell 中运行 `ipconfig`，找到 **vEthernet (WSL)** 适配器的 IPv4 地址，通常以 `172.` 开头，例如 `172.31.0.1`。

> 如果 WSL 网络模式为 **Mirrored**（WSL Settings → 网络 → 网络模式），则宿主机 IP 可以直接使用 `127.0.0.1`（即 localhost），无需查找 vEthernet 地址。

### 3. 修改 Clash DNS 监听

检查 Clash 配置文件，将 DNS 监听地址从 `127.0.0.1:53` 改为 `0.0.0.0:53`，允许 WSL 通过宿主机 IP 访问 DNS。

---

## 三、WSL 中配置代理

在 `~/.bashrc` 末尾添加 `proxy_on` 和 `proxy_off` 函数，方便随时开关代理：

```bash
proxy_on() {
  HOST_IP="127.0.0.1"
  export http_proxy="http://$HOST_IP:7890"
  export https_proxy="http://$HOST_IP:7890"
  git config --global http.proxy "http://$HOST_IP:7890"
  git config --global https.proxy "http://$HOST_IP:7890"
  echo "Proxy ON: $HOST_IP:7890"
}

proxy_off() {
  unset http_proxy https_proxy
  git config --global --unset http.proxy
  git config --global --unset https.proxy
  echo "Proxy OFF"
}
```

保存后 `source ~/.bashrc` 生效。使用时：

```bash
proxy_on    # 开启代理（安装软件、git clone 前）
proxy_off   # 关闭代理（日常使用、连接集群时）
```

> **注意**：`HOST_IP` 根据你的 WSL 网络模式填写。Mirrored 模式用 `127.0.0.1`，NAT 模式用 vEthernet 的 `172.x.x.x` 地址。

---

## 四、VSCode 安装插件及连接 WSL

### 1. 推荐插件

在 VSCode 扩展商店（`Ctrl+Shift+X`）搜索并安装：

**基础必备：**

| 插件 | 用途 |
|---|---|
| **Chinese (Simplified) Language Pack** (`ms-ceintl.vscode-language-pack-zh-hans`) | 中文界面汉化 |
| **Remote - WSL** (`ms-vscode-remote.remote-wsl`) | 连接 WSL 子系统 |
| **Remote - SSH** (`ms-vscode-remote.remote-ssh`) | 远程连接集群/服务器 |

**Git 与协作：**

| 插件 | 用途 |
|---|---|
| **GitLens** (`eamodio.gitlens`) | Git 增强：行级 blame、历史对比、代码作者追踪 |
| **Git Graph** (`mhutchie.git-graph`) | Git 提交历史可视化 |

**AI 辅助：**

| 插件 | 用途 |
|---|---|
| **Claude Code** (Anthropic 官方扩展) | AI 辅助编程 |

**编辑增强：**

| 插件 | 用途 |
|---|---|
| **Markdown All in One** (`yzhang.markdown-all-in-one`) | Markdown 编辑增强（快捷键、目录生成、预览） |
| **Markdown Preview Mermaid Support** (`bierner.markdown-mermaid`) | 在 Markdown 预览中渲染 Mermaid 图表 |
| **Path Intellisense** (`christian-kohler.path-intellisense`) | 文件路径自动补全 |
| **Todo Tree** (`gruntfuggly.todo-tree`) | 汇总并高亮代码中的 TODO/FIXME 标记 |

### 2. 连接 WSL

三种方式任选：

- **方式一**：WSL 终端中进入项目目录，输入 `code .`，VSCode 会自动以 Remote-WSL 模式打开
- **方式二**：VSCode 左下角点击绿色 `><` 图标，选择 **"连接到 WSL"**
- **方式三**：`Ctrl+Shift+P` → 输入 `WSL` → 选择 **"Remote-WSL: New Window"**

连接成功后左下角显示 **"WSL: Ubuntu"**，终端自动切换为 WSL 的 bash。

---

## 五、WSL 中安装开发工具

> 以下操作均在 WSL 终端中执行。如果安装过程中卡住，先执行 `proxy_on` 开启代理。

### 1. 安装 nvm 和 Node.js

```bash
# 安装 nvm
curl -o- https://raw.githubusercontent.com/nvm-sh/nvm/v0.39.7/install.sh | bash
source ~/.bashrc

# 安装 Node.js LTS 版本
nvm install 22
nvm use 22
nvm alias default 22

# 验证
node --version   # v22.x.x
npm --version    # 10.x.x
```

### 2. 安装 Claude Code

```bash
npm install -g @anthropic-ai/claude-code
```

安装完成后 `claude` 命令即可在终端中使用。同时建议在 VSCode 扩展商店安装 Claude Code 插件，获得更好的编辑体验。

### 3. 安装 ccswitch

[ccswitch](https://github.com/TomokiMatsubuchi/ccswitch) 是一个 Claude Code 配置切换工具，用于管理不同的 API 提供商配置。

从 GitHub Releases 页面下载 **CC-Switch-v3.16.2-Windows.msi**，双击安装。安装完成后在 Windows 开始菜单中即可找到 CC-Switch。

---

## 六、ccswitch 接入 DeepSeek

### 1. 新增 DeepSeek 供应商

在 Windows 中启动 CC-Switch，在界面中新增一个 DeepSeek 供应商，填写以下信息：

| 配置项 | 值 |
|---|---|
| API Base URL | `https://api.deepseek.com/anthropic` |
| API Key | `sk-你的DeepSeek_API_Key` |
| 默认模型 | `deepseek-v4-pro`（或 `deepseek-v4-flash`） |

CC-Switch 会将此配置写入 Claude Code 的配置目录。

### 2. 修改 VSCode 中 Claude Code 的设置

在 VSCode 中打开 Claude Code 扩展的设置：

- **设置 → 通用 → 窗口行为**：根据需要调整（如"在新标签页打开"）
- **设置 → 高级 → 配置文件目录**：指向 ccswitch 生成的配置目录（通常为 `~/.claude/`）

这两个设置确保 VSCode 中的 Claude Code 扩展使用与终端中 `claude` 命令相同的 DeepSeek 后端配置。

### 3. 验证配置

配置完成后，`~/.claude/settings.json` 应包含类似内容：

```json
{
  "env": {
    "ANTHROPIC_AUTH_TOKEN": "sk-xxx",
    "ANTHROPIC_BASE_URL": "https://api.deepseek.com/anthropic",
    "ANTHROPIC_MODEL": "deepseek-v4-pro",
    "ANTHROPIC_DEFAULT_HAIKU_MODEL": "deepseek-v4-flash",
    "ANTHROPIC_DEFAULT_OPUS_MODEL": "deepseek-v4-pro",
    "ANTHROPIC_DEFAULT_SONNET_MODEL": "deepseek-v4-pro"
  }
}
```

> DeepSeek 的 Anthropic 兼容接口将 Haiku 映射为 `deepseek-v4-flash`（快速模型），Sonnet/Opus 映射为 `deepseek-v4-pro`（推理模型）。可根据实际需要调整。

---

## 七、VSCode 中使用 Claude Code

### 1. 打开 Claude Code

- **终端模式**：WSL 终端中直接输入 `claude`，进入交互式对话
- **侧边栏模式**：VSCode 左侧点击 Claude Code 图标，或 `Ctrl+Shift+P` → **"Claude Code: Open"**

### 2. 基本用法

| 操作 | 说明 |
|---|---|
| 直接提问 | 在对话框中输入问题，Claude 会读取当前项目文件进行回答 |
| 代码编辑 | Claude 可直接修改文件，修改前会展示 diff 供确认 |
| 终端命令 | Claude 可执行 bash 命令，需用户批准 |
| `@file` 引用 | 使用 `@文件名` 让 Claude 聚焦特定文件 |
| 选中代码提问 | 在编辑器中选中代码段，`Ctrl+Shift+L` 发送到 Claude |

### 3. 常用 Slash 命令

| 命令 | 作用 |
|---|---|
| `/help` | 查看帮助 |
| `/fast` | 切换快速模式（使用更快/更便宜的模型） |
| `/config` | 修改配置（主题、模型等） |
| `/clear` | 清空对话历史 |
| `/memory` | 管理持久化记忆 |

---

## 八、Claude Code 安装 Skills

Skills 是 Claude Code 的功能扩展包，可以为特定领域（科学计算、文档处理、可视化等）提供专用能力。

### 1. 查看已安装的 Skills

```bash
ls ~/.claude/skills/
```

### 2. 安装 Skills

- **通过 marketplace**：在 Claude Code 对话中输入 `/plugin` 打开 marketplace，浏览并安装。K-Dense AI Claude Scientific Skills 等科学计算系列 skills 均可在 marketplace 中找到
- **通过 skill-creator**：使用 `/skill-creator` 自定义创建专用 skill
- **通过 find-skills**：使用 `/find-skills 关键词` 按需搜索，例如 `/find-skills 分子动力学`

### 3. 推荐 Skills

**系统与开发：**

| Skill | 用途 |
|---|---|
| `superpowers` | 通用能力增强：规划、执行、测试、调试一站式工作流 |
| `using-superpowers` | superpowers 的使用指南 |
| `skill-creator` | 自定义创建专用 skill |
| `find-skills` | 按关键词搜索可用 skills |
| `dispatching-parallel-agents` | 并行处理多个独立任务 |
| `subagent-driven-development` | 子代理驱动的开发流程 |
| `using-git-worktrees` | Git 多分支管理 |

**科学计算（K-Dense AI Claude Scientific Skills 系列）：**

| Skill | 用途 |
|---|---|
| `pymatgen` / `ase` | 材料科学建模与分析 |
| `scientific-writing` | 学术写作辅助 |
| `scientific-visualization` | 科学数据可视化 |
| `scientific-brainstorming` | 科研头脑风暴 |
| `scientific-critical-thinking` | 科研批判性思维 |
| `scientific-schematics` | 科学原理图/示意图绘制 |
| `scientific-slides` | 学术演示文稿制作 |
| `literature-review` | 文献综述辅助 |
| `hypothesis-generation` | 科研假设生成 |
| `peer-review` | 模拟同行评审 |
| `molecular-dynamics` | 分子动力学模拟 |
| `matplotlib` / `seaborn` | 数据绘图 |
| `polars` / `pandas-like` | 高性能数据处理 |
| `rdkit` / `datamol` / `molfeat` | 化学信息学与分子特征 |
| `transformers` / `pytorch-lightning` | 深度学习 |

**文档与办公：**

| Skill | 用途 |
|---|---|
| `xlsx` / `pptx` / `pdf` / `docx` | 办公文档处理 |
| `markdown-mermaid-writing` | Markdown 中绘制 Mermaid 图表 |
| `infographics` | 信息图制作 |
| `latex-posters` | LaTeX 学术海报 |

---

## 九、Skills 使用小技巧

1. **不用记名字**：输入 `/find-skills` 按关键词搜索，例如 `/find-skills 绘图` 会列出所有可视化相关 skills
2. **技能即命令**：安装后直接 `/skill-name` 调用，skill 会自动加载对应领域的上下文和工具
3. **按需加载**：skill 在首次调用时才加载，不会拖慢启动速度，不必担心装太多
4. **组合使用**：例如先用 `xlsx` 处理数据，再用 `matplotlib` 绘图，最后用 `pptx` 插入到演示文稿
5. **自定义 skill**：如果有重复性的工作流（如固定格式的 VASP 输入文件生成），用 `/skill-creator` 创建一个专用 skill，后续一键调用
6. **更新 skill**：通过 marketplace 检查更新，保持 skill 与上游同步

---

## 十、快速排错

| 问题 | 原因 | 解决方法 |
|---|---|---|
| `npm install` 卡住 | 未开代理 | 执行 `proxy_on` 后重试 |
| `claude` 命令不存在 | npm 全局路径未加入 PATH | 确认 `~/.nvm/versions/node/v22.22.0/bin` 在 PATH 中 |
| Claude Code 无法连接 API | 代理拦截或 API Key 无效 | `proxy_off` 关闭代理后重试；检查 `~/.claude/settings.json` 中的 `ANTHROPIC_AUTH_TOKEN` |
| VSCode 无法连接 WSL | WSL 未启动 | PowerShell 中执行 `wsl` 启动子系统 |
| ccswitch 切换后不生效 | VSCode 扩展未重启 | `Ctrl+Shift+P` → "Reload Window" |
