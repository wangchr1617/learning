# 晶格热导率计算

可参考教程：1. [B站视频](https://www.bilibili.com/video/BV1wo4y1U772?t=30.4) 2. [Phonopy官网](https://phonopy.github.io/phonopy/index.html) 3. ShengBTE+thirdorder手册

## 1. 高精度优化（VASP）

为保证优化收敛，可以精度逐渐调高，先同时优化晶格和原子位置用`ISIF = 3`，后面只优化原子位置用`ISIF = 2`。逐渐增大`EDIFFG`的的精度。
INCAR参数设置：

    EDIFFG = -1E-8
    ISIF = 2
    IBRION = 1

完美晶格优化具体步骤可以看wb/thermalconductivity/perfect文件夹下的opt2-opt5（opt1就是当前文件夹下的vasp计算文件）文件夹中的`INCAR`。（金刚石没那么难优化，其他的复杂体系可以按照这样的优化步骤。）

## 2. 有限位移法计算获得二阶力常数（Phonopy+VASP）

##### 1. 对优化过后的晶胞扩胞

    $ phonopy -d --dim='2 2 2'

##### 2. 计算二阶力常数

对生成的POSCAR-\*文件全部进行静态计算，建立band.conf。

    $ phonopy -f vasprun.xml-00*
    $ phonopy --dim='2 2 2' -c POSCAR-unitcell -p -s band.conf
    $ phonopy-bandplot --gnuplot>222.dat

生成`FORCE_CONSTANTS`
完美晶格`FORCE_CONSTANTS`文件夹位置:
wb/thermalconductivity/perfect/FC\_2ND/FORCE\_CONSTANTS
空位晶格同理。

计算声子谱，

1.  可以使用IBRION=5/6/7/8,5和6是使用有限位移法进行计算，7和8是使用紧束缚近似方法（DFPT）方法进行计算，紧束缚无法使用VDW。
2.  使用任何以上一个参数计算完成后，建立声子谱输入文件band.conf

<!---->

    ATOM_NAME = Ge Se 
    DIM = 5 5 1
    BAND = 0.0 0.0 0.0  0.333 0.333 0.0  0.5 0.0 0.0  0.0 0.0 0.0
    BAND_LABELS = G K M G
    FORCE_CONSTANTS = READ
    FC_SYMMETRY = .TRUE.

生成力常数文件

    phonopy --fc vasprun.xml

    phonopy -f 1/vasprun.xml 2/vasprun.xml 3/vasprun.xml 4/vasprun.xml 5/vasprun.xml 6/vasprun.xml 7/vasprun.xml 8/vasprun.xml 9/vasprun.xml 10/vasprun.xml 11/vasprun.xml 12/vasprun.xml 13/vasprun.xml 14/vasprun.xml 15/vasprun.xml 16/vasprun.xml 17/vasprun.xml 18/vasprun.xml 19/vasprun.xml 20/vasprun.xml 21/vasprun.xml 22/vasprun.xml 23/vasprun.xml 24/vasprun.xml 25/vasprun.xml 26/vasprun.xml 27/vasprun.xml 28/vasprun.xml 29/vasprun.xml 30/vasprun.xml

1.  使用该命令直接输出声子谱pdf图片

<!---->

    phonopy -c POSCAR-unitcell -p -s band.conf

声子谱详计算细信息：<https://zhuanlan.zhihu.com/p/599295736>

**ps**：修正，band.conf中，设置了`FC_SYMMETRY=.TRUE.`考虑了对称性，因此生成的`FORCE_CONSTANT`第一行两个数字（移动位置的原子数，扩胞原子数）并不相等。需将第二行命令更改为：

    $ phonopy --dim='2 2 2' -c POSCAR-unitcell band.conf --full-fc

生成不含对称性的`FORCE_CONSTANT`，不然使用ShengBTE的时候会报错。

## 3. 三阶力常数矩阵计算（thirdorder+VASP）

##### 1. 对优化过后的晶胞扩胞

    $ thirdorder_vasp.py sow 2 2 2 -8

sow a b c m  abc表示生成的超胞大小，m为正数的时候表示截断半径（单位：nm）（例如0.6，表示截断半径为6.0埃）；为负数时表示第几近邻截断（此时应取负整数，例如-3表示第三近邻）对所有生成的POSCAR都进行高精度静态计算，INCAR与有限位移法计算二阶力常数保持一致。

##### 2. 计算三阶力常数

    $ find job* -name vasprun.xml|sort -n|thirdorder_vasp.py reap 2 2 2 -8

生成`FORCE_CONSTANTS_3RD`
完美晶格`FORCE_CONSTANTS_3RD`文件位置：
wb/thermalconductivity/perfect/FC\_3RD/FORCE\_CONSTANTS\_3RD
空位晶格同理。

**ps**：chatgpt4.0对于力常数的解释：
三阶力常数是指在固体物理学和材料科学中用来描述原子间相互作用随它们相对位置的三次变化而变化的参数。这是一个高阶的概念，相对于更常见的二阶力常数而言，它涉及到更复杂的相互作用。

在晶格动力学中，二阶力常数涉及原子间作用力随两个原子之间的位移的二次变化，而三阶力常数则涉及作用力随位移的三次变化。具体来说：

1.  二阶力常数：描述了当一个原子从其平衡位置微小位移时，作用于它的力如何随着这个位移线性变化（即力与位移成正比）。

2.  三阶力常数：在此基础上，考虑原子位移更复杂的效应，即当原子从平衡位置发生位移时，作用于它的力不仅与位移成正比，还与位移的平方成比例。这意味着原子间的相互作用不仅取决于它们之间的距离，还取决于这些距离变化的方式。

三阶力常数对于理解和描述非线性效应至关重要，比如非线性晶体振动、热导率的非线性行为以及材料的非线性声学性质。在热力学性质的计算、非线性声子相互作用的分析以及高压下物质行为的研究中，三阶力常数是一个重要的参数。通过这些高阶力常数，科学家能更准确地预测和解释复杂材料系统中的物理现象。

## 4. 晶格热导率计算（ShengBTE）

1.  将`CONTROL`，`FORCE_CONSTANTS`，`FORCE_CONSTANTS_3RD`放到一个文件夹下，将`FORCE_CONSTANTS`改为`FORCE_CONSTANTS_2ND`。
2.  提交ShengBTE.pbs进行计算。
3.  一共计算了100K-700K（间隔100K）七个温度下的晶格热导率。晶格热导率在每个温度文件夹中的TiK（i代表温度）文件夹下的`BTE.kappa_tensor`文件的最后一行。每列对应的方向顺序分别为XX,XY,XZ,YX,YY,YZ,ZX,ZY,ZZ。取XX，YY，ZZ作为X，Y，Z方向的晶格热导率值。例如完美晶格在300K下的晶格热导率的文件位置在wb/thermalconductivity/perfect/300K/T300K/BTE.kappa\_tensor。空位同理。
4.  所有软件中的文件信息可以查阅软件说明书文件夹中的 ShengBTE-README.md（无效链接）。

# 电子热导率计算（BoltzTrap2+VASP）

热导率*kappae* = 晶格热导率*kappae*(l)+电子热导率*kappae*(e)

##### 计算电子热导率步骤：

1.  结构优化，静态计算（使用的是和计算晶格热导率一样的K点分布）
2.  在静态计算文件夹中使用BoltzTrap2:

    \$ btp2 -vv interpolate -m 3 wb/perfect/scf/case

生成interpolation.bt2文件

    $ btp2 -vv integrate interpolation.bt2 200:700:50

得到的interpolation.trace记录了体系输运因子的信息。

1.  interpolation.trace中第8列`kappae/tau0[W/(m*K*s)]`是可以用来计算电子热导率的数据

计算电子热导率的公式，以及公式中需要用到的参数

# 热容计算（Phonopy）

参考教程：1. [Phonopy官网](https://phonopy.github.io/phonopy/index.html) 2. [知乎教程](https://zhuanlan.zhihu.com/p/448370499)

1.  在计算二阶力常数的文件夹中建立mesh.conf文件

2.  计算热力学性质

    \$ phonopy -t mesh.conf -c POSCAR-unitcell

3.  生成热力学dat文件方便作图

    \$ phonopy-propplot --gnuplot thermal\_properties.yaml > thermal.dat

完美晶格文件在wb/thermalconductivity/perfect/FC\_2ND/thermal.dat。空位晶格同理。

# 热膨胀系数（Phonopy+VASP）

计算步骤：

1.  结构优化得到CONTCAR改成POSCAR

2.  改写缩放系数0.95-1.05建立11个POSCAR

3.  对11个POSCAR进行声子计算

4.  计算热学性质并收集thermal\_properties-{1..11}.yaml文件

5.  提取e-v.dat文件

6.  phonopy-qha计算热学性质。

完美晶格：
步骤见wb/thermalexpansivity/perfect/thermalexpansion.sh
热膨胀系数dat文件：wb/thermalexpansivity/perfect/CP/thermal\_expansion.dat
空位晶格同理。

# 杨氏模量，泊松比计算（vaspkit+VASP）

教程：[公众号文章](https://mp.weixin.qq.com/s?__biz=MzI2OTQ4OTExOA==\&mid=2247487112\&idx=1\&sn=bfa1e8c7981b15e880cd8426c7f8a8ab\&chksm=eadec839dda9412f12c9020b189fd767a1d9d1bfbc1929822c9410bc8c24834a330b663b0038\&mpshare=1\&scene=1\&srcid=1104qk2jJQfIvkhhJR8u318F\&sharer_sharetime=1604472889425\&sharer_shareid=dbb0734284422c0aebeb989a7ec537c0#rd)

1.  构建原胞。
2.  高精度优化。
3.  准备VPKIT.in文件，使用vaspkit 200生成发生应变的结构文件。
4.  使用步骤2中的INCAR对新生成的结构进行计算。
5.  将VPKIT.in文件第一行改为2，使用vaspkit 200计算。
6.  可以将弹性常数矩阵放入[ELATE](https://progs.coudert.name/elate)分析。

力学数据保存在：
完美晶格：wb/Youngmodulus/perfect/opt/opt2/data.txt
空位晶格：wb/Youngmodulus/vacancy1/opt/data.txt

# 光学性质计算（vaspkit+VASP）

教程：
[知乎帖子1](https://www.zhihu.com/question/287951798)
[知乎帖子2](https://zhuanlan.zhihu.com/p/669635613)
[公众号文章](https://mp.weixin.qq.com/s/or3KbvwAv_h3B68AJAB1fA)

#### 方法：

1.  结构优化
2.  静态计算得到WAVECAR和CHGCAR
3.  使用静态计算的WAVECAR和CHGCAR进行光学计算
    注意NBANDS的取值一般为自洽计算中OUTCAR的NBANDS的值的2-3倍。
4.  运行vaspkit 711得到光学性质

完美晶格结果在wb/optics/perfect/dielectricconstant
光吸收系数：ABSORPTION.dat
能量损失谱：ENERGY\_LOSSSPECTRUM.dat
消光系数：EXTINCTION.dat
复介电函数虚部：IMAG.in（无效链接）
复介电函数实部：REAL.in（无效链接）
反射系数：REFLECTIVITY.dat
折射系数：REFRACTIVE.dat

空位晶格同理。
输出文件 ABSORB.dat，REFRACTIVE.dat，REFLECTIVITY.dat，EXTINCTION.dat 和 ENERGYLOSSSPECTRUM.dat，依次为 absorption coefficient, refractive coefficient, reflectivity coefficient, extinction coefficient and energy-loss function。
**ps**：

1.  第一次计算后发现光学性质的曲线毛刺非常多，可能是3.中的`NEDOS`的值没有取高，建立test2文件夹设置`NEDOS = 2000`进行测试。
2.  精度仍然较差，建立test3文件夹设置`NEDOS = 10000`进行测试。
3.  增大2. 中的K点密度。生成新的WAVECAR和CHGCAR。建立test4文件夹进行测试。曲线变得比较平滑。所以光学性质参考test4中的数据。

