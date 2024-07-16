import numpy as np
import matplotlib.pyplot as plt
from gpyumd.load import load_dos, load_vac

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

num_corr_steps = 1000
dos = load_dos(num_dos_points=num_corr_steps)['run0']
vac = load_vac(num_corr_steps)['run0']
dos['DOSxyz'] = dos['DOSx'] + dos['DOSy'] + dos['DOSz']
vac['VACxyz'] = vac['VACx'] + vac['VACy'] + vac['VACz']
vac['VACxyz'] /= vac['VACxyz'].max()
print('DOS:', list(dos.keys()))
print('VAC:', list(vac.keys()))

plt.figure(figsize=(12, 10))
plt.subplot(2, 2, 1)
plt.plot(vac['t'], vac['VACx'] / vac['VACx'].max(), color='C3', linewidth=3)
plt.plot(vac['t'], vac['VACy'] / vac['VACy'].max(), color='C0', linestyle='--', linewidth=3)
plt.plot(vac['t'], vac['VACz'] / vac['VACz'].max(), color='C2', linestyle='-.', zorder=100, linewidth=3)
plt.xlabel('Correlation Time (ps)')
plt.ylabel('VAC (Normalized)')
plt.xlim(0, 5)
plt.legend(['x', 'y', 'z'])
plt.title('(a)')

plt.subplot(2, 2, 2)
plt.plot(dos['nu'], dos['DOSx'], color='C3', linewidth=3)
plt.plot(dos['nu'], dos['DOSy'], color='C0', linestyle='--', linewidth=3)
plt.plot(dos['nu'], dos['DOSz'], color='C2', linestyle='-.', zorder=100, linewidth=3)
plt.xlabel(r'$\nu$ (THz)')
plt.ylabel('PDOS (1/THz)')
plt.xlim(0, 8)
plt.ylim(0, )
plt.legend(['x','y', 'z'])
plt.title('(b)')

plt.subplot(2, 2, 3)
plt.plot(vac['t'], vac['VACxyz'], color='k', linewidth=3)
plt.xlabel('Correlation Time (ps)')
plt.ylabel('VAC (Normalized)')
plt.xlim(0, 5)
plt.title('(c)')

plt.subplot(2, 2, 4)
plt.plot(dos['nu'], dos['DOSxyz'], color='k', linewidth=3)
plt.xlabel(r'$\nu$ (THz)')
plt.ylabel('PDOS (1/THz)')
plt.xlim(0, 8)
plt.ylim(0, )
plt.title('(d)')

plt.tight_layout()
plt.savefig('./VDOS.png', bbox_inches='tight')
