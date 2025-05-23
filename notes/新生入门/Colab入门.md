
# Colab 入门

Colaboratory 是一个免费的 Jupyter 笔记本环境，不需要进行任何设置就可以使用，并且完全在云端运行。
借助 Colaboratory，可以编写和执行代码、保存和共享分析结果，以及利用强大的计算资源，所有这些都可通过浏览器免费使用。

## 如何使用 Colab

Colab 一般配合 Google Drive 使用，访问 [谷歌云端硬盘](https://workspace.google.com/products/drive/#download) 下载 GoogleDriveSetup.exe 并安装 Google Drive。
打开安装好的 Google Drive 应用，打开 `我的云端硬盘`：

<div align="left">
<img src="./figures/colab_001.png" width = "50%" />
</div>

在 `新建` → `更多` 里单机右键创建 `Google Colaboratory`，如图所示：

<div align="left">
<img src="./figures/colab_002.png" width = "50%" />
</div>

如果没有 `Google Colaboratory` 选项，那就选择 `关联更多应用` 并搜索 `Colaboratory` 下载。

新建的 `.ipynb` 文件如图所示：

<div align="left">
<img src="./figures/colab_003.png" width = "50%" />
</div>

可以像 Jupyter Notebook 一样使用 Colab，甚至他们的文件格式都是一样的。
Jupyter Notebook 文件可以直接传到谷歌云并用 Colab 打开并执行：

<div align="left">
<img src="./figures/colab_004.png" width = "50%" />
</div>

如果要使用 GPU 资源，需要选择 `Runtime` 并 `Change runtime type`，如图所示：

<div align="left">
<img src="./figures/colab_005.png" width = "50%" />
</div>

<div align="left">
<img src="./figures/colab_006.png" width = "50%" />
</div>

使用以下代码即可查看当前获取的 GPU 资源：
```
from google.colab import drive 
drive.mount('/content/drive') 

gpu_info = !nvidia-smi 
gpu_info = '\n'.join(gpu_info) 
if gpu_info.find('failed') >= 0: 
    print('Not connected to a GPU') 
else: 
    print(gpu_info) 
```
代码执行情况如图所示：
<div align="left">
<img src="./figures/colab_007.png" width = "50%" />
</div>