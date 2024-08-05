"""
用法: python vasp_xdatcar_analysis.py

功能描述:
该脚本用于处理VASP的XDATCAR输出文件，计算声子态密度（Phonon DOS）、速度自相关函数（VACF），
以及对一对原子类型的径向分布函数（Pair Correlation Function, PCF）。

输出文件:
- pdos.dat: 声子态密度
- vacf.dat: 速度自相关函数
- pcf.dat: 径向分布函数

输出图像:
- pdos.png: 声子态密度图
- vacf.png: 速度自相关函数图
- pcf.png: 径向分布函数图

注意事项:
请确保在INCAR文件中设置NBLOCK=1，以便所有配置都写入XDATCAR。
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.fftpack import fft, fftfreq
from scipy.ndimage import gaussian_filter1d

# 常量定义
kb = 8.617332478E-5  # Boltzmann Constant in [eV/K]
ev = 1.60217733E-19  # electron volt in [Joule]
Navogadro = 6.0221412927E23  # Avogadro's Constant

class XDATCAR:
    """Python Class for VASP XDATCAR file analysis."""

    def __init__(self, file='XDATCAR'):
        """
        初始化XDATCAR类。

        参数:
        file (str): XDATCAR文件路径，默认值为'XDATCAR'。
        """
        self.xdatcar = file

        # 初始化属性
        self.potim = None  # MD时间步长
        self.mtype = None  # 每种原子类型的质量
        self.read_outcar()  # 从OUTCAR中读取POTIM和POMASS

        self.TypeName = None
        self.ChemSymb = None
        self.Ntype = None
        self.Nions = None
        self.Nelem = None
        self.Niter = None

        # Direct Coordinate中的位置
        self.position = None
        # Cartesian Coordinate中的位置
        self.positionC = None
        # 速度，单位为 Å/fs
        self.velocity = None
        self.read_xdat()

        self.mass_and_name_per_ion()
        # 温度
        self.Temp = np.zeros(self.Niter - 1)
        # 动能
        self.Ken = np.zeros(self.Niter - 1)
        # 时间，单位为fs
        self.Time = np.arange(self.Niter - 1) * self.potim
        self.get_temp()

        # 速度自相关函数
        self.VAF = None
        self.VAF2 = None

    def mass_and_name_per_ion(self):
        """计算每个离子的质量和化学符号。"""
        self.mions = []
        self.ChemSymb = []

        if self.TypeName is None:
            # 使用字母A-Z表示原子类型名称
            self.TypeName = [chr(i) for i in range(65, 91)][:self.Ntype]

        for i in range(self.Ntype):
            self.mions += [np.tile(self.mtype[i], self.Nelem[i])]
            self.ChemSymb += [np.tile(self.TypeName[i], self.Nelem[i])]

        self.mions = np.concatenate(self.mions)
        self.ChemSymb = np.concatenate(self.ChemSymb)

    def read_xdat(self):
        """读取VASP XDATCAR文件以获取位置和速度数据。"""
        inp = [line for line in open(self.xdatcar) if line.strip()]
        scale = float(inp[1])
        
        self.cell = np.array([line.split() for line in inp[2:5]], dtype=float)
        self.cell *= scale

        ta = inp[5].split()
        tb = inp[6].split()
        
        if ta[0].isalpha():
            self.TypeName = ta
            self.Ntype = len(ta)
            self.Nelem = np.array(tb, dtype=int)
            self.Nions = self.Nelem.sum()
        else:
            # 处理VASP 4.X格式
            self.Nelem = np.array(tb, dtype=int)
            self.Nions = self.Nelem.sum()
            self.Ntype = len(tb)
            self.TypeName = None

        pos = np.array([line.split() for line in inp[7:] if not line.split()[0].isalpha()], dtype=float)
        self.position = pos.ravel().reshape((-1, self.Nions, 3))
        self.Niter = self.position.shape[0]

        dpos = np.diff(self.position, axis=0)
        self.positionC = np.zeros_like(self.position)
        
        # 应用周期性边界条件
        dpos[dpos > 0.5] -= 1.0
        dpos[dpos < -0.5] += 1.0
        
        # 计算速度，单位为 Å/fs
        for i in range(self.Niter - 1):
            self.positionC[i, :, :] = np.dot(self.position[i, :, :], self.cell)
            dpos[i, :, :] = np.dot(dpos[i, :, :], self.cell) / self.potim

        self.positionC[-1, :, :] = np.dot(self.position[-1, :, :], self.cell)
        self.velocity = dpos

    def read_outcar(self):
        """从OUTCAR文件中读取POTIM和POMASS信息。"""
        if os.path.isfile("OUTCAR"):
            # print "OUTCAR found!"
            # print "Reading POTIM & POMASS from OUTCAR..."

            outcar = [line.strip() for line in open('OUTCAR')]
            lp = 0
            lm = 0
            for ll, line in enumerate(outcar):
                if 'POTIM' in line:
                    lp = ll
                if 'Mass of Ions in am' in line:
                    lm = ll + 1
                if lp and lm:
                    break

            # 从OUTCAR中获取POTIM和每种原子类型的质量
            self.potim = float(outcar[lp].split()[2])
            self.mtype = np.array(outcar[lm].split()[2:], dtype=float)

    def get_temp(self, Nfree=None):
        """计算温度随时间的变化。

        参数:
        Nfree (int): 自由度，默认值为 3 * (Nions - 1)。
        """
        for i in range(self.Niter - 1):
            ke = np.sum(np.sum(self.velocity[i, :, :] ** 2, axis=1) * self.mions / 2.0)
            self.Ken[i] = ke * 1E7 / Navogadro / ev
            if Nfree is None:
                Nfree = 3 * (self.Nions - 1)
            self.Temp[i] = 2 * self.Ken[i] / (kb * Nfree)

    def get_vaf(self):
        """计算速度自相关函数（VAF）。"""
        self.VAF2 = np.zeros((self.Niter - 1) * 2 - 1)
        
        for i in range(self.Nions):
            for j in range(3):
                self.VAF2 += np.correlate(self.velocity[:, i, j], self.velocity[:, i, j], 'full')
        
        # 两侧VAF
        self.VAF2 /= np.sum(self.velocity ** 2)
        self.VAF = self.VAF2[self.Niter - 2:]

    def phonon_dos(self, unit='THz', sigma=5):
        """从VAF计算声子态密度（Phonon DOS）。

        参数:
        unit (str): 频率单位，默认值为 'THz'。
        sigma (int): 高斯平滑参数，默认值为5。

        返回:
        tuple: 包含频率和态密度的元组。
        """
        N = self.Niter - 1
        
        # 频率转换
        omega = fftfreq(2 * N - 1, self.potim) * 1E3
        
        if unit.lower() == 'cm-1':
            omega *= 33.35640951981521
        elif unit.lower() == 'mev':
            omega *= 4.13567
        
        if self.VAF2 is None:
            self.get_vaf()
        
        # 高斯平滑
        smVAF = gaussian_filter1d(self.VAF2, sigma=sigma)
        pdos = np.abs(fft(smVAF - np.average(smVAF))) ** 2

        return omega[:N], pdos[:N]

    def pair_correlation_function(self, bins=50, Niter=10, A='', B=''):
        """计算径向分布函数（Pair Correlation Function, PCF）。

        参数:
        bins (int): 直方图的bins数量，默认值为50。
        Niter (int): 迭代步数，默认值为10。
        A (str): 第一种元素类型，默认值为第一个类型。
        B (str): 第二种元素类型，默认值为第一个类型。

        返回:
        tuple: 包含径向分布函数值和bins边界的元组。
        """
        if not A:
            A = self.TypeName[0]
        if not B:
            B = A

        whichA = self.ChemSymb == A
        whichB = self.ChemSymb == B
        indexA = np.arange(self.Nions)[whichA]
        indexB = np.arange(self.Nions)[whichB]
        posA = self.position[:, whichA, :]
        posB = self.position[:, whichB, :]

        steps = range(0, self.Niter, Niter)
        rABs = np.array([
            posA[i, k, :] - posB[i, j, :]
            for k in range(indexA.size)
            for j in range(indexB.size)
            for i in steps
            if indexA[k] != indexB[j]
        ])
        
        # 应用周期性边界条件
        rABs[rABs > 0.5] -= 1.0
        rABs[rABs < -0.5] += 1.0

        # 从直坐标转换到笛卡尔坐标
        rABs = np.linalg.norm(np.dot(self.cell, rABs.T), axis=0)

        # 计算对间距的直方图
        val, b = np.histogram(rABs, bins=bins)
        
        # 计算系统的密度
        rho = self.Nions / np.linalg.det(self.cell)
        # A类型原子的数量
        Na = self.Nelem[self.TypeName.index(A)]
        # B类型原子的数量
        Nb = self.Nelem[self.TypeName.index(B)]
        dr = b[1] - b[0]
        
        # 归一化处理
        val = val * self.Nions / (4 * np.pi * b[1:] ** 2 * dr) / (Na * Nb * rho) / len(steps)

        return val, b[1:]

if __name__ == '__main__':
    # 创建XDATCAR对象并计算相关函数
    inp = XDATCAR('./XDATCAR')
    inp.get_vaf()
    
    # 计算声子态密度
    unit = 'cm-1'    
    omega, dos = inp.phonon_dos(unit=unit)
    np.savetxt('pdos.dat', np.array([omega, dos]).T)
    
    # 绘制声子态密度图
    plt.plot(omega, dos)
    plt.xlabel('Frequency (' + unit + ')')
    plt.ylabel('DOS (u)')
    plt.xlim((0, 5000))
    plt.savefig('pdos.png', dpi=1000)
    plt.close() 
    plt.clf()

    # 计算并保存速度自相关函数
    np.savetxt('vacf.dat', np.array([inp.Time, inp.VAF]).T)
    plt.plot(inp.Time, inp.VAF)
    plt.ylabel('VACF')
    plt.xlabel('Time (fs)')
    plt.xlim((0, 1000))
    plt.savefig('vacf.png', dpi=1000)
    plt.close()
    plt.clf()

    # 计算并保存径向分布函数
    val, bins = inp.pair_correlation_function(bins=1000, A='O', B='H')
    np.savetxt('pcf.dat', np.array([bins, val]).T)
    plt.plot(bins, val)
    plt.ylabel('g(r)')
    plt.xlabel('r (Å)')
    plt.savefig('pcf.png', dpi=1000)
