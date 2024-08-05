"""
用法: python tdos.py

功能描述:
该脚本读取 TDOS.dat 文件中的数据，并生成总态密度（TDOS）图。
结果图像保存为 TDOS.png 文件。

注意事项:
请确保 TDOS.dat 文件位于脚本运行的当前目录中。
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# 读取 TDOS.dat 文件
with open('./TDOS.dat', 'r') as f:
    data = []
    for line in f.readlines():
        # 将每一行拆分为列表并过滤掉空字符串
        items = [item for item in line.replace("\n", "").split(" ") if item != ""]
        data.append(items)

# 将数据转换为 Pandas DataFrame
df = pd.DataFrame(data)

# 转换列为数值类型，忽略第一行（假设第一行是表头）
x = pd.to_numeric(df[0][1:].reset_index(drop=True), errors='coerce')
y = pd.to_numeric(df[1][1:].reset_index(drop=True), errors='coerce')

# 绘制总态密度 (TDOS) 图
plt.plot(x, y, color="r", linewidth=1, label="TDOS")
plt.xlabel('Energy (eV)')
plt.ylabel('DOS')
plt.xlim(x.min(), x.max())
plt.ylim(0, y.max() * 1.1)
plt.legend(loc='best')

# 保存图像为 TDOS.png，设置bbox_inches='tight'以去除多余边距
plt.savefig('./TDOS.png', bbox_inches='tight')
