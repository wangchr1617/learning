
# VASP 弹性常数计算 - 有限位移法

INCAR参考设置如下：

```
ISTART = 0
ICHARG = 2
LCHARG = .F.
LWAVE = .F.

ALGO = Normal
NELM = 120
EDIFF = 1E-08		# 有限差分法要求提高收敛精度

ISIF = 3			# 弹性常数计算时 ISIF >= 3
IBRION = 6			# 这里 5 表示振动频率计算，6 表示弹性常数计算
NSW = 1				# 只进行一步离子弛豫
NFREE = 4			# 这里 2 和 4 都可以，它们的差距很小
POTIM = 0.015		# 在 IBRION = 5 / 6 时，POTIM 决定离子移动的步长，单位为埃

ENCUT = 800			# 计算弹性常数时一般需要较高的平面波截断能
ISMEAR = 0
SIGMA = 0.05
PREC = Accurate
LREAL = .False.
ADDGRID = .TRUE.
```

由于计算弹性常数时对 K 点密度要求较高，使用 `vaspkit 102` 生成 KPOINTS 时 K-mesh 应该小于等于 0.03。
注意，不要设置 `NPAR`、`NCORE` 等并行参数，否则会报错。

计算完成后，弹性常数矩阵会直接写入 OUTCAR。
可以使用命令 `grep -A9 ‘TOTAL ELASTIC MODULI’ OUTCAR` 查看弹性刚度矩阵（注意单位 kBar = 0.1 GPa）。
对不同的晶系的晶体，因为对称性的关系，弹性刚度矩阵中独立的弹性常数是确定的。
因此，晶系的对称性越高，独立的张量元数目越少。

利用 `vaspkit 203` 可以直接得到弹性常数和机械性能。
得到弹性常数后，体积弹性模量剪切模量和杨氏模量的估计值是 Voigh-Reuss-Hill（VRH）的平均值，维氏硬度是用 Chen-Niu 模型计算得到的，而断裂韧性是基于 Niu-Niu-Oganov 模型计算得到的。
