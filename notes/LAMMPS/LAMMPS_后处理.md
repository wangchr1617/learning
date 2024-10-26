
# LAMMPS 后处理

## 提取热力学量

LAMMPS 的 in 文件中，使用 `thermo_modify line yaml` 将自定义的热力学量输出为 YAML 格式，
以方便 python 从日志文件 `log.lammps` 中提取和解析这些数据：
```
import re, yaml
import pandas as pd
import matplotlib.pyplot as plt

plt.rcParams['font.size'] = 14
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
	
try:
    from yaml import CSafeLoader as Loader
except ImportError:
    from yaml import SafeLoader as Loader
 
docs = ""
with open("log.lammps") as f:
    for line in f:
        m = re.search(r"^(keywords:.*$|data:$|---$|\.\.\.$|  - \[.*\]$)", line)
        if m: docs += m.group(0) + '\n'
 
thermo = list(yaml.load_all(docs, Loader=Loader))
 
print("Number of runs: ", len(thermo))
print("Columns of run 0: ", thermo[0]['keywords'][:])

print("-----------------------------")

df = pd.DataFrame(data=thermo[1]['data'], columns=thermo[1]['keywords'])
fig = df.plot(x='Temp', y=['Cella', 'Cellb', 'Cellc',], ylabel='Cell length')
plt.savefig('cell_length.png', bbox_inches='tight')
```