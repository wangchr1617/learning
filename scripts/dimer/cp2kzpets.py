#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
此脚本通过读入做频率计算的cp2k.out文件
输出频率的值，单位：cm-1
输出-TS的值，单位：eV
输出ZPE的值，单位：eV
输出ZPE - TS的值，单位：eV
"""

import sys
import math
from scipy import constants as con
import numpy as np

import getopt

opts,args = getopt.getopt(sys.argv[1:],'-h',['help'])
for opt_name,opt_value in opts:
    if opt_name in ('-h','--help'):
        print('''此脚本通过读入做频率计算的cp2k.out文件
输出频率的值，单位：cm-1
输出-TS的值，单位：eV
输出ZPE的值，单位：eV
输出ZPE - TS的值，单位：eV''')
        exit()



f = open('cp2k.out','r')
freq0 = []
lines = f.readlines()
for line in lines:
   if "VIB|Frequency" in line:
     freq0.append(float(line.strip('\n').split()[2]))
     freq0.append(float(line.strip('\n').split()[3]))
     freq0.append(float(line.strip('\n').split()[4]))
for i in range(len(freq0)):
    if 0 < freq0[i] < 50:
        freq0[i] = 50
#num_ignore是需要去掉的频率
#目的是只留下振动频率
#对表面吸附分子：0，线性分子(如双原子分子O2)：5，非线性分子：6
script, num_ignore = sys.argv
#num_ignore = "0"
num_ignore = int(num_ignore)
#freq0 = [420.609746, 402.026105]
freq1 = np.array(freq0) #读取频率计算中各个频率的值，单位：cm-1
freq1 = freq1[num_ignore:] #忽略num_ignore个频率，留下振动频率
freq2 = freq1[freq1 > 0] #计算自由能时需要把虚频去除

ZPE = (freq2.sum()) * 0.123984 / 2000 #eV



h_p = con.h   # 6.62606957E-34 # J*s Plank Constant
k_b = con.k   # 1.38064852E-23 # m2kg*s-2*k-1 Boltzman Constant
R_gas = con.R # 8.3144598      # J*mol-1*K-1 Gas Constant
l_s = con.c   # 299792458      # light speed m * s-1
Tem = 298.15     # Temperature    # K
beta = 1/(k_b * Tem)

def get_pf(nu): # get partition function
    x_i = h_p * float(nu) * l_s * beta
    pf_l = x_i / (math.exp(x_i) - 1) # Left part in the entropy equation
    pf_r = math.log(1 - math.exp(-x_i))
    pf   = pf_l - pf_r
    entropy = R_gas * pf
    return entropy
#script, nu = sys.argv
#nu = 402.026105
entropy = []
ts = []
for nu in freq2:
    nu = float(nu) * 100 #Convert cm-1 to m-1
     
    entropy0 = get_pf(nu) # J * K-1 * mol-1
    ts0      = entropy0 * Tem / 1000 / 96.485 # in eV
    entropy.append(entropy0)
    ts.append(ts0)
entropy_sum = np.array(entropy).sum()
ts_sum = np.array(ts).sum()
print("Frequency (cm-1):")
for fre in freq0:
    print(fre)
print(("Total number of frequencies: " + str(len(freq0))))
print(("The number of frequencies used when calculating ZPE-TS: " + str(len(freq2))))
print("********************************")
print("-TS (eV):")
print((-ts_sum))

print("ZPE (eV)")
print(ZPE)
print(" ")
print("---**************************---")
print("ZPE - TS (eV)")
print((ZPE - ts_sum))
