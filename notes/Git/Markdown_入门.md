
# Markdown的使用说明

Markdown 是一种轻量级的标记语言，能够将纯文本格式转换为 HTML 页面。它的语法简单明了，易于阅读和编写。

### 标题
使用 `#` 来表示标题，共有六级标题，可以根据需要选择使用：

```markdown
# 一级标题
## 二级标题
### 三级标题
#### 四级标题
##### 五级标题
###### 六级标题
```

### 列表
有序列表和无序列表的使用：

```markdown
- 无序列表项1
- 无序列表项2
    - 嵌套无序列表项

1. 有序列表项1
2. 有序列表项2
    1. 嵌套有序列表项
```

### 强调
使用 `*` 或 `_` 表示斜体，使用 `**` 或 `__` 表示加粗：

```markdown
*斜体*
_斜体_
**加粗**
__加粗__
```

### 插入链接

使用 `[描述](链接)` 的格式来插入链接：

```markdown
[Markdown官方文档](https://www.markdownguide.org)
```

### 插入图片

使用 `![描述](图片链接)` 的格式来插入图片：

```markdown
![Markdown Logo](https://markdown-here.com/img/icon256.png)
```
或者使用 html 语言：
```html
<div align="left">
<img src="./figures/fig_001.png" width = "50%" />
</div>
```

### 插入表格

使用 `|` 来创建表格，并使用 `-` 来分隔表头和表格内容：

```markdown
| 头1 | 头2 | 头3 |
| --- | --- | --- |
| 内容1 | 内容2 | 内容3 |
| 内容4 | 内容5 | 内容6 |
```
或者使用 html 语言：
```html
<table>
    <tr>
        <th>姓名</th>
        <th>数学</th>
        <th>英语</th>
        <th>科学</th>
    </tr>
    <tr>
        <td>张三</td>
        <td>90</td>
        <td>85</td>
        <td>88</td>
    </tr>
    <tr>
        <td>李四</td>
        <td>78</td>
        <td>90</td>
        <td>85</td>
    </tr>
    <tr>
        <td>王五</td>
        <td>85</td>
        <td>80</td>
        <td>92</td>
    </tr>
</table>
```

使用 Markdown 语法编辑表格时，可以使用一些在线工具快速对齐表格，例如 [Markdown Tables generator](https://www.tablesgenerator.com/markdown_tables) 等。

### 插入文件

Markdown 不支持直接插入文件，但可以通过链接的方式来提供文件的下载：

```markdown
[下载文件](链接地址)
```

### 插入代码

可以使用反引号 `` ` `` 来插入行内代码，使用三个反引号来插入代码块：

```markdown
行内代码示例：`print("Hello, World!")`

代码块示例：
\`\`\`python
def hello_world():
    print("Hello, World!")

hello_world()
\`\`\`
```

### 插入 LaTeX 数学公式

Markdown 支持插入 LaTeX 数学公式。行内公式使用 `$`，块级公式使用 `$$`：

```markdown
行内公式示例：$E = mc^2$

块级公式示例：
$$
E = mc^2
$$
```

1. **常用符号**

| 名称  | 示例                 | 显示                       |
|-----|--------------------|--------------------------|
| 上下标 | A=x_1^2+x_2^2      | $A=x_1^2+x_2^2$          |
| 根号  | \sqrt{3} 或 \sqrt[3]{x} | $\sqrt{3}$ 或 $\sqrt[3]{x}$ |    
| 分数  | \frac{1}{2} 或 \dfrac{x}{y} | $\frac{1}{2}$ 或 $\dfrac{x}{y}$ |            

2. **特殊运算符号**

| markdown代码 | 显示        | markdown代码 | 显示 |
|------------|-----------|-|-|
| \pm        | $\pm$     | \mp | $\mp$ |
| \div       | $\div$    | \times | $\times$ |
| \geq       | $\geq$    | \leq | $\leq$ |
| \approx    | $\approx$ | \neq | $\neq$ |
| \sum       | $\sum$    | \equiv | $\equiv$ |
| \log       | $\log$    | \prod | $\prod$ |
| \lg        | $\lg$     | \ln | $\ln$ |
| \int       | $\int$    | \sin | $\sin$ |
| \lim       | $\lim$    | \iint | $\iint$ |
| \notin     | $\notin$  | \in | $\in$ |
| \subset    | $\subset$ | | |

3. **希腊字母**

| markdown代码 | 显示         | markdown代码 | 显示        |
|------------|------------|------------|-----------|
| \alpha     | $\alpha$   | \beta      | $\beta$   |
| \Gamma     | $\Gamma$   | \gamma     | $\gamma$  |
| \Delta     | $\Delta$   | \delta     | $\delta$  |
| \epsilon   | $\epsilon$ | \zeta      | $\zeta$   |
| \eta       | $\eta$     | \Theta     | $\Theta$  |
| \theta     | $\theta$   | \Lambda    | $\Lambda$ |
| \lambda    | $\lambda$  | \mu        | $\mu$     |
| \nu        | $\nu$      | \xi        | $\xi$     |
| \pi        | $\pi$      | \rho       | $\rho$    |
| \Sigma     | $\Sigma$   | \sigma     | $\sigma$  |
| \Phi       | $\Phi$     | \phi       | $\phi$    |
| \chi       | $\chi$     | \Psi       | $\Psi$    |
| \psi       | $\psi$     | \Omega     | $\Omega$  |
| \omega     | $\omega$   | \tau       | $\tau$    |
