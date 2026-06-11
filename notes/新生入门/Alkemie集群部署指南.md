# Alkemie 远程集群部署 + 本地浏览器访问指南

## 环境信息

| 项目 | 值 |
|---|---|
| 集群地址 | `changruiwang-ICME@10.251.0.28`，端口 2022 |
| 集群系统 | Linux（mu01 节点） |
| 集群内存 | 62GB 总 / 55GB 可用 |
| Python 环境 | conda 环境 `alkemie`，Python 3.12.12 |
| 本地机器 | macOS（Apple Silicon） |
| 服务端口 | 8445（8443 被占用时更换） |

---

## 一、本地构建与打包

### 1. 构建前端

```bash
cd /path/to/Alkemie
npm install
npm run build
mv build backend/build
```

### 2. 打包代码

```bash
tar --exclude='node_modules' \
    --exclude='.git' \
    --exclude='__pycache__' \
    --exclude='*.pyc' \
    --exclude='venv' \
    --exclude='.venv' \
    -czf alkemie-deploy.tar.gz \
    backend/ public/ src/ package.json
```

### 3. 上传到集群

```bash
scp -P 2022 alkemie-deploy.tar.gz changruiwang-ICME@10.251.0.28:~/
```

---

## 二、集群上部署

### 1. 解压并安装依赖

```bash
ssh -p 2022 changruiwang-ICME@10.251.0.28
mkdir -p ~/Software/Alkemie && cd ~/Software/Alkemie
tar xzf ~/alkemie-deploy.tar.gz
conda activate alkemie

# 分批安装，避免峰值内存叠加（--no-cache-dir 减少磁盘/内存占用）
pip install --no-cache-dir flask flask-cors flask-socketio python-socketio \
    eventlet python-dotenv werkzeug pyyaml markdown pydantic httpx

pip install --no-cache-dir numpy
pip install --no-cache-dir scipy
pip install --no-cache-dir ase
pip install --no-cache-dir pymatgen
pip install --no-cache-dir phonopy seekpath
pip install --no-cache-dir openai sentence-transformers
pip install --no-cache-dir sqlalchemy pymongo
pip install --no-cache-dir psycopg2-binary PyJWT bcrypt pytest
pip install --no-cache-dir dscribe scikit-learn umap-learn dpdata
```

### 2. 配置 .env

```bash
cp backend/.env.example backend/.env
vim backend/.env
```

关键配置：
```ini
ALKEMIE_WORKSPACES_ROOT=~/alkemie-workspaces
ALKEMIE_ALLOWED_ORIGINS=http://localhost:8445,http://localhost:3000
LLM_PROVIDER=deepseek
LLM_API_KEY=sk-你的key
LLM_BASE_URL=https://api.deepseek.com
LLM_MODEL=deepseek-v4-flash
HF_HUB_OFFLINE=1
TRANSFORMERS_OFFLINE=1
```

> `ALKEMIE_ALLOWED_ORIGINS` 必须包含 `http://localhost:8445`，否则浏览器刷新后会报 `not an accepted origin` 错误。

### 3. 启动服务

```bash
cd ~/Software/Alkemie/backend
conda activate alkemie
ALKEMIE_PORT=8445 python app.py
```

成功日志：
```
INFO:__main__:Background monitor started
INFO:websocket_manager:Starting WebSocket server on 0.0.0.0:8445
```

---

## 三、本地浏览器访问

### 1. 建立 SSH 隧道

```bash
ssh -L 8445:localhost:8445 -p 2022 changruiwang-ICME@10.251.0.28 -N
```

> `-N` 表示不执行命令，只维持隧道。终端卡住不动 = 隧道正常工作。

### 2. 浏览器访问

打开 `http://localhost:8445`

> 注意：是 `http://` 不是 `https://`。集群未配置 SSL 证书，使用纯 HTTP。

---

## 四、遇到的困难与解决方法

| 序号 | 问题 | 原因 | 解决方法 |
|---|---|---|---|
| 1 | `Address already in use` (8443) | 端口被之前残留的进程占用 | `fuser -k 8443/tcp` 或换端口 |
| 2 | `Killed` (进程被杀死) | OOM Killer 内存超限 | 可能是旧进程 OOM 导致 socket 未释放；加 `PYTHONMALLOC=malloc` 降低内存分配器开销；分批安装依赖避免峰值 |
| 3 | `Address already in use` (8444) | 旧进程被 kill 后 socket 仍处于 TIME_WAIT | 换全新端口 8445 |
| 4 | 浏览器 `https://` 无法打开 | 集群未配置 SSL 证书，服务是纯 HTTP | 改用 `http://localhost:8445`|
| 5 | SSH 隧道命令卡住 | `-N` 参数的正确行为，隧道在运行中 | 不需要操作，保持终端开着即可 |
| 6 | WebSocket disconnected 页面刷新后出现 | 集群进程挂掉或 SSH 隧道断开 | 先 `fuser -k 8445/tcp` 杀残留进程，重启服务，确认隧道存活 |
| 7 | `not an accepted origin` (WebSocket 连接被拒) | `.env` 中 `ALKEMIE_ALLOWED_ORIGINS` 未包含 `http://localhost:8445` | 编辑 `backend/.env`，确保 `ALKEMIE_ALLOWED_ORIGINS=http://localhost:8445` |
| 8 | `Ctrl+C` 无法停止服务 | eventlet 协程未响应 SIGINT | 另开终端 `fuser -k 8445/tcp` 强制杀进程 |

---

## 五、安全退出

### 1. 停止集群上的服务

在运行 `python app.py` 的终端中按 `Ctrl+C` 即可正常退出。

> 如果 `Ctrl+C` 无效（eventlet 协程有时不响应 SIGINT），另开一个 SSH 终端强制杀进程：
> ```bash
> ssh -p 2022 changruiwang-ICME@10.251.0.28 "fuser -k 8445/tcp"
> ```

不要直接关终端，否则进程可能残留在后台，下次启动时端口被占用。

如果终端已丢失，需要手动杀进程：

```bash
ssh -p 2022 changruiwang-ICME@10.251.0.28
fuser -k 8445/tcp
```

### 2. 断开本地 SSH 隧道

在运行 SSH 隧道的终端中按 `Ctrl+C` 即可断开。

### 3. 关闭浏览器

直接关闭标签页即可，无持久状态需要保存。

---

## 六、快速重启（跳过部署）

部署后的服务文件已持久化在集群磁盘上，日常使用或演示无需重新构建、上传或安装。每次启动只需两步：

### 第一步：启动集群服务

```bash
ssh -p 2022 changruiwang-ICME@10.251.0.28
cd ~/Software/Alkemie/backend && conda activate alkemie
ALKEMIE_PORT=8445 python app.py
```

看到以下日志表示启动成功：
```
INFO:websocket_manager:Starting WebSocket server on 0.0.0.0:8445
```

### 第二步：本地建立隧道 + 浏览器访问

本地新开终端：

```bash
ssh -L 8445:localhost:8445 -p 2022 changruiwang-ICME@10.251.0.28 -N
```

浏览器打开 `http://localhost:8445`。

> 全程约 10-15 秒（取决于集群上 Python 导入耗时）。

---

## 七、更新代码后重新部署

当 wchr 分支有代码更新时需要重新部署：

```bash
# 本地重新构建前端 + 打包
cd /path/to/Alkemie
npm run build && mv build backend/build
tar --exclude='node_modules' --exclude='.git' --exclude='__pycache__' \
    --exclude='*.pyc' --exclude='venv' --exclude='.venv' \
    -czf alkemie-deploy.tar.gz backend/
scp -P 2022 alkemie-deploy.tar.gz changruiwang-ICME@10.251.0.28:~/

# 集群上：先停服务（Ctrl+C），再解压覆盖
tar xzf ~/alkemie-deploy.tar.gz -C ~/Software/Alkemie/
# 然后按第六章步骤重新启动
```
