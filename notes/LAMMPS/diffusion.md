#### 仅输出空位周围的原子

首先计算出每个原子的配位数（指定截断）；然后通过dump命令，输出那些配位数小于完美晶格的原子。

```
compute coord all coord/atom cutoff 2.8
dump events all custom 1 dump.prd id type x y z
dump_modify events thresh c_coord != 4
```

#### 仅输出发生了空位扩散时的扩散原子

```
write_dump      all custom tmp.dump id type x y z    # 先输出一次所有原子的结构，后续输出的第一帧位移为0，是空的。

variable        Dhop equal 0.6
variable        check atom "c_dsp[4] > v_Dhop"
compute         dsp all displace/atom refresh check
dump            1 all custom 100 tmp.dump id type x y z
dump_modify     1 append yes thresh c_dsp[4] > ${Dhop} refresh c_dsp delay 100
```

