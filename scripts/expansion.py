import matplotlib.pyplot as plt
import numpy as np
from gpyumd.load import load_thermo  # 假设此为正确导入方式

aw = 2
fs = 16
font = {'size': fs}
plt.rc('font', **font)
plt.rc('axes', linewidth=aw)
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'

def set_fig_properties(ax_list):
    tl = 8
    tw = 2
    tlm = 4

    for ax in ax_list:
        ax.tick_params(which='major', length=tl, width=tw)
        ax.tick_params(which='minor', length=tlm, width=tw)
        ax.tick_params(which='both', axis='both', direction='in', right=True, top=True)

thermo = load_thermo()
print("Thermo quantities:", list(thermo.keys()))

t = 0.01 * np.arange(1, thermo['temperature'].shape[0] + 1) 
Natoms = 12096  # Number of atoms
NT = 1000  # dump_thermo
temp = np.array([300, 400, 500, 550, 560, 570,
                 580, 590, 600, 610, 620, 630,
                 640, 650, 700, 800])
M = thermo['temperature'].shape[0] // NT
v_ave = (thermo['Lx'] * thermo['Ly'] * thermo['Lz']) / Natoms
Pave = (thermo['Px'] + thermo['Py'] + thermo['Pz']) / 3.

plt.figure(figsize=(12, 10))
plt.subplot(2, 2, 1)
plt.plot(t, thermo['temperature'])
plt.xlim(0, t.max())  
plt.ylim(thermo['temperature'].min()-50, thermo['temperature'].max()+50) 
plt.xlabel('Time (ps)')
plt.ylabel('Temperature (K)')
plt.title('(a)')

plt.subplot(2, 2, 2)
plt.plot(t, Pave)
plt.xlim(0, t.max()) 
plt.ylim(-0.2, 0.2)  
plt.yticks(np.arange(-0.2, 0.25, 0.1))  
plt.xlabel('Time (ps)')
plt.ylabel('Pressure (GPa)')
plt.title('(b)')

plt.subplot(2, 2, 3)
plt.plot(t, v_ave, linewidth=3)
plt.xlim(0, t.max())  
plt.ylim(v_ave.min()-0.1, v_ave.max()+0.1)  
plt.xlabel('Time (ps)')
plt.ylabel(r'Volume (Å³/atom)')
plt.yticks(np.arange(int(v_ave.min()), int(v_ave.max())+1, 0.2))  
plt.title('(c)')

plt.subplot(2, 2, 4)
plt.scatter(temp, v_ave[::100], s=20, zorder=100, facecolor='none', edgecolors='C0', linewidths=3)
def fit(x,y):
    coefficients = np.polyfit(x, y, 1)
    slope = coefficients[0]
    intercept = coefficients[1]
    plt.plot(x, slope*x + intercept, 'r-')
fit(temp[:4], v_ave[::100][:4])
fit(temp[-4:], v_ave[::100][-4:])
plt.xlim(thermo['temperature'].min()-50, thermo['temperature'].max()+50) 
plt.ylim(v_ave.min()-0.1, v_ave.max()+0.1)  
plt.xlabel('Temperature (K)')
plt.ylabel(r'Volume (Å³/atom)')
plt.yticks(np.arange(int(v_ave.min()), int(v_ave.max())+1, 0.2))  
plt.title('(d)')

plt.tight_layout() 
plt.savefig('./expansion.png', bbox_inches='tight')

