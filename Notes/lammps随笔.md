
# lammps 随笔

---

### 后处理技巧


LAMMPS 支持热样式“yaml”，对于“自定义”样式热力学输出，可以使用thermo_modify line yaml将格式更改为YAML。

使用 python 能方便地从日志文件中提取和解析这些数据：
```
import re, yaml
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
print(thermo[1]['keywords'][4], ' = ', thermo[1]['data'][2][4])
```