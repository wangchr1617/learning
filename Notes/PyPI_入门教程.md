# PyPI 入门

---

### 注册账号

- 在 [测试 PyPI](https://test.pypi.org/) 平台注册账号，用于测试发布的库。
- 在 [正式 PyPI](https://pypi.org/) 平台注册账号，用于正式发布库。

记得申请 API 接口，并获取 token。

---

### 创建项目

创建一个新的目录来作为库的项目文件夹。这个目录中需要包含一个 `setup.py` 文件来描述你的库，并且在适当的位置创建 Python 模块文件以及其他必要的文件。

```
altbc_analyzer/             # 目录
├── altbc_analyzer          # 子目录
│   ├── altbc_analyzer.py   # 自定义 Python 类一
│   ├── __init__.py         # 将当前目录标记为 Python 模块
│   └── neighbor_list.py    # 自定义 Python 类二
├── analyze_and_plot.py
├── example
│   └── POSCAR
├── LICENSE                 # 许可证
├── README.md               # 说明文档
├── requirements.txt        # 依赖库
└── setup.py                # 安装脚本
```

---

### 编写代码

为了后续 `import` 方便，我们在 `__init__.py` 里对外部暴露的包名规范一下，例如：

```python
from .altbc_analyzer import ALTBC_Analyzer
from .neighbor_list import NeighborList
```

---

### 编写 README.md 和 LICENSE

- `README.md` 遵循 markdown 语法，介绍项目的基本信息、安装方法、使用示例等。
- `LICENSE` 在 GitHub 创建项目时即可选择，确保项目的版权信息。

---

### 编写 setup.py

首先编写发布库的 `setup.py` 文件。`setup.py` 是 `setuptools` 的构建脚本，用于描述项目。它告诉 PyPI 项目名称、版本、依赖库、支持的操作系统、支持的 Python 版本等基本信息。手写 `setup.py` 难度较大，个人建议 `git clone  https://github.com/kennethreitz/setup.py` ，使用仓库里的 `setup.py` 模板，并修改必要的配置即可。

下面是 `altbc_analyzer` 项目中我修改的部分：

```python
# Package meta-data.
NAME = 'altbc_analyzer'
DESCRIPTION = 'A package for calculating ALTBC in crystal structures.'
URL = 'https://github.com/wangchr1617/altbc_analyzer'
EMAIL = 'wangchr1617@gmail.com'
AUTHOR = 'Changrui Wang'
REQUIRES_PYTHON = '>=3.6.0'
VERSION = '0.1.0'

# What packages are required for this module to be executed?
REQUIRED = [
    'ase',
    'numpy',
    'pandas',
    'matplotlib',
    'scipy',
]
```

---

### 发布到 PyPI 测试环境

一旦确认库已经完善，可以选择将其发布。为了保证 `pip install` 正常，建议先把库发布到 PyPI 的测试环境。

#### Windows

Windows 系统使用 PowerShell 脚本 `upload_pypi_test.ps1` 上传到 PyPI 的测试环境：

```powershell
Remove-Item -Recurse -Force .uild
Remove-Item -Recurse -Force .\dist
Remove-Item -Recurse -Force .ltbc_analyzer.egg-info

python setup.py sdist bdist_wheel

$env:TWINE_USERNAME = "__token__"
$env:TWINE_PASSWORD = "XXX XXX XXX"

python -m twine upload --repository-url https://test.pypi.org/legacy/ dist/*
```

注意这里使用了前面申请的 API token。

如果执行脚本遇到权限问题，可能需要调整 PowerShell 的执行策略。可以在 PowerShell 中运行以下命令来允许运行脚本：

```powershell
Set-ExecutionPolicy RemoteSigned
```

或者对于本地脚本，可以使用：

```powershell
Set-ExecutionPolicy Bypass -Scope Process
```

#### Linux

在 Linux 系统中，改用下面这个 Shell 脚本 `upload_pypi_test.sh` 上传即可：

```sh
rm -rf ./build
rm -rf ./dist
rm -rf ./altbc_analyzer.egg-info

python setup.py sdist bdist_wheel

export TWINE_USERNAME="__token__"
export TWINE_PASSWORD="XXX XXX XXX"

python -m twine upload --repository-url https://test.pypi.org/legacy/ dist/*
```

发布到 PyPI 测试环境后，可以直接通过 `pip install -i https://test.pypi.org/simple/ altbc_analyzer` 命令安装。在 Python 中成功导入即可准备正式发布。

---

### 正式发布到 PyPI 平台

只需将上面的上传脚本中的

```sh
python -m twine upload --repository-url https://test.pypi.org/legacy/ dist/*
```

替换为

```sh
python -m twine upload dist/*
```

即可用于正式发布。
