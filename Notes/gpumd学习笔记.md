# GPUMD & NEP

2024-07-16

---

### 樊老师讲义

1. 经验势

2. NEP 机器学习势

3. 经典分子动力学模拟

4. 路径积分分子动力学

---

### 前、后处理仓库安装
Windows 系统中安装 calorine 时，注意修改`calorine/src/nepy/nepy.cpp`文件中的
```
#include <unistd.h>
```
为
```
#ifdef _WIN32
#include <windows.h>
#else
#include <unistd.h>
#endif
```
即可`pip install .`成功安装 calorine。

---

### 训练集生成策略

参见 [GPUMD-wizard.ipynb]

### 单点计算和训练集提取

---

### 训练集筛选策略

1. PCA + FPS

2. UMAP + FPS

---

### 主动学习

1. 基于不确定性的主动学习

- MCMD 探索

- 温度探索

- 单轴拉伸

---

### NEP 训练

1. NEP 训练可视化

2. NEP + ZBL

3. NEP + DFTD3

4. 偶极矩和电极化张量NEP

---

### NEP 的基准测试

1. Wizard 性质验证

2. RDF & ADF

---

### GPUMD 模拟

1. 弛豫

2. 拉伸

3. 辐照损伤

---

### 基于 GPUMD 的性质模拟

1. 自扩散系数

2. 粘滞系数

3. 热导率

4. 热导率校正

5. 谱热导率

6. 声子态密度 (SHC)

7. 声子平均自由程 (MFP)

8. κALDo + NEP

9. Dynasor轨迹分析

10. ISAACS轨迹分析

---