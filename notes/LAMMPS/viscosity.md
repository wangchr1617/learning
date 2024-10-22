## G-K法计算粘度

### 学习lammps示例（二维）

```
# settings

variable        x equal 20
variable        y equal 20

variable        rho equal 0.6							# 晶格常数
variable        t equal 1.0								# 温度
variable        rc equal 2.5							# 势函数截断半径

variable    p equal 400     							# correlation length 关联长度
variable    s equal 5       							# sample interval 取样间隔
variable    d equal $p*$s   							# dump interval dump文件储存间隔

# problem setup

units           lj
dimension       2
atom_style      atomic
neigh_modify    delay 0 every 1

lattice         sq2 ${rho}
region          simbox block 0 $x 0 $y -0.1 0.1
create_box      1 simbox
create_atoms    1 box

pair_style      lj/cut ${rc}
pair_coeff      * * 1 1

mass            * 1.0
velocity        all create $t 97287

# equilibration run

fix             1 all nve
fix             2 all langevin $t $t 0.1 498094
fix             3 all enforce2d							# 二维材料

thermo          $d
run             10000

velocity        all scale $t

unfix           2

# Green-Kubo viscosity calculation

reset_timestep  0

# Define distinct components of symmetric traceless stress tensor

variable         pxy equal pxy
variable         pxx equal pxx-press

# ave/correlate在不同的时间间隔计算它们之间的时间相关性，并在较长的时间尺度上平均相关数据；并记录到文件profile.gk.2d中 利用ave running 计算出的数据是平均值
fix              SS all ave/correlate $s $p $d v_pxy v_pxx type auto file profile.gk.2d ave running

# Diagonal components of SS are larger by factor 2-2/d,
# which is 4/3 for d=3, but 1 for d=2.
# See Daivis and Evans, J.Chem.Phys, 100, 541-547 (1994)

variable         scale equal 1.0/$t*vol*$s*dt
variable         diagfac equal 2-2/2
variable         vxy equal trap(f_SS[3])*${scale}
variable         vxx equal trap(f_SS[4])*${scale}/${diagfac}

thermo_style     custom step temp press pxy v_vxy v_vxx

run              500000

variable         etaxy equal v_vxy
variable         etaxx equal v_vxx
variable         eta equal 0.5*(${etaxy}+${etaxx})
print            "running average viscosity: ${eta}"
```

## 对于三维流体

