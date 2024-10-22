## 1. 裂纹扩展知识

在真实材料中裂纹扩展速率非常快，和声速相当！

## 2. 模拟结构

![image-20231012181306255](C:\Users\zefengli\AppData\Roaming\Typora\typora-user-images\image-20231012181306255.png)

## 3. lammps输入文件

```lammps
# Crack simulation
dimension 3
units metal
boundary s s p
atom_style atomic

# Create geometry
lattice fcc 3.52
region box block 0 100 0 30 -5 5
create_box 5 box						# 5 指有五种类型的原子
create_atoms 1 box						# 用 1 类型原子填满盒子

# Set mass for atom types
mass 1 58.69
mass 2 58.69
mass 3 58.69
mass 4 58.69
mass 5 58.69

# EAM potentials
pair_style eam
pair_coeff * * Ni_u3.eam

# Define groups
region 1 block INF INF INF 1.25 INF INF					# 还是用的lattice单位
group lower region 1
region 2 block INF INF 28.75 INF INF INF
group upper region 2
group boundary union lower upper						# 将lower和upper两个组的原子组合为一个组，然后用总的减去边界就是中间的可移动原子
group mobile subtract all boundary

region leftupper block INF 20 15 INF INF INF
region leftlower block INF 20 INF 15 INF INF
group leftupper region leftupper
group leftlower region leftlower

# Set atom types
set group leftupper type 2
set group leftlower type 3
set group lower type 4
set group upper type 5

# Compute stress and Voronoi volume
compute s all stress/atom NULL					# 计算每个原子单位体积内的应力，分为6个方向。
compute vol all voronoi/atom
variable von_press atom sqrt(0.5*((c_s[1]-c_s[2])^2+(c_s[1]-c_s[3])^2+(c_s[2]-c_s[3])^2+6*(c_s[4]^2+c_s[5]^2+c_s[6]^2)))/10000/c_vol[1]

# Variables for strain and stress				# 计算Y方向应力应变的变化
variable l equal ly								# ly 内置的变量，l会随着ly变化
variable len equal ${l}							# Y方向初始长度，值不会变化
variable strain equal (ly-v_len)/v_len
variable stress equal -pyy/10000 #convert from bar to GPa		# pyy 是内置变量，y方向应力，单位从bar转换为GPa

# Initial velocities
compute new mobile temp							# 只统计mobile区域的温度
velocity mobile create 1 887723 temp new		# mobile区域速度初始化
velocity upper set 0.0 0.1 0.0					# upper组向上速度，0.1埃/ps
velocity mobile ramp vy 0 0.1 y 1.25 28.75 sum yes				# 沿y方向线性变化

# Fixes
fix 1 all nve									# 为什么用NVE？
fix 2 boundary setforce NULL 0.0 0.0			# 对X方向受力不做干预

# Run simulation
timestep 0.001
thermo 1000								# 每隔1000步输出热力学信息
thermo_modify temp new					# 输出中的温度只统计mobile部分
neigh_modify exclude type 2 3			# 实现裂纹萌生。紧邻列表排除23，则这两部分原子间没有相互作用

# Output
dump 1 mobile custom 1000 dump.crack id type x y z v_von_press
fix 3 all ave/time 5 40 1000 v_strain v_stress  file ori.crack		# 每隔1000步，取最后的200步每隔5步采样一次，一共进行40次采样，然后把采样得到的平均应力应变输出到指定文件
run 100000
```

### 3.1 残余应力计算

$$
\bar{\sigma}=\sqrt{\frac{1}{2}\left[\left(\sigma_x-\sigma_y\right)^2+\left(\sigma_x-\sigma_z\right)^2+\left(\sigma_y-\sigma_z\right)^2+6\left(\sigma_{x y}^2+\sigma_{x z}^2+\sigma_{y z}^2\right)\right]}
$$

其中σ_x/σ_y/σz为主应力，σxy/σyz/σxz为剪切应力。

*注意这里计算残余应力是每个原子上的残余应力。* 

```lammps
# Compute stress and Voronoi volume
compute s all stress/atom NULL					# 计算每个原子单位体积内的应力，分为6个方向。该方法就是将原子附近的一个立方体作为单元，计算其力；实际上计算出来的是应力*单位体积：P*V
compute vol all voronoi/atom					# 计算每个原子的体积
variable von_press atom sqrt(0.5*((c_s[1]-c_s[2])^2+(c_s[1]-c_s[3])^2+(c_s[2]-c_s[3])^2+6*(c_s[4]^2+c_s[5]^2+c_s[6]^2)))/10000/c_vol[1]		# 10000指bar到GPa的单位换算
```

## 4. 问题

1. 在24核的batch节点运行时，报错：ERROR on proc 2: Divide by 0 in variable formula。但是在32核的manycores节点运行没有报错。为什么？

这是脚本中采用了 **Voronoi/atom** 方法计算原子体积，该方法可能会将某个原子的体积置为0，导致程序崩溃。

对于存在大量真空的体系，采用 **Voronoi/atom** 计算原子体积也会有问题，因为所有原子的体积加起来相当于盒子的体积。这相当于把真空层也用来计算体积，计算出的原子体积是不准确的。

## 5. 后续分析

### 5.1 位错密度统计

