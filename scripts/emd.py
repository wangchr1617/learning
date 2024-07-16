import numpy as np
import matplotlib.pyplot as plt
from gpyumd.load import load_hac

hac = load_hac([50000]*3, [10]*3)
print("Runs:", list(hac.keys()))
print("Run Data:", list(hac['run0'].keys()))

t = hac['run0']['t']
hac_ave_i = np.zeros(hac['run0']['jxijx'].shape[0])
hac_ave_o = np.zeros_like(hac_ave_i)
ki_ave, ko_ave = np.zeros_like(hac_ave_i), np.zeros_like(hac_ave_o)
for runkey in hac.keys():
    hac_ave_i += hac[runkey]['jxijx'] + hac[runkey]['jyijy']
    hac_ave_o += hac[runkey]['jxojx'] + hac[runkey]['jyojy']
    ki_ave += (hac[runkey]['kxi'] + hac[runkey]['kyi'])
    ko_ave += (hac[runkey]['kxo'] + hac[runkey]['kyo'])
hac_ave_i /= hac_ave_i.max()
hac_ave_o /= hac_ave_o.max()
ki_ave /= 6.
ko_ave /= 6.

aw = 2
fs = 16
font = {'size': fs}
plt.rc('font', **font)
plt.rc('axes', linewidth=aw)

def set_fig_properties(ax_list):
    tl = 8
    tw = 2
    tlm = 4

    for ax in ax_list:
        ax.tick_params(which='major', length=tl, width=tw)
        ax.tick_params(which='minor', length=tlm, width=tw)
        ax.tick_params(which='both', axis='both', direction='in', right=True, top=True)

plt.figure(figsize=(12, 4))

plt.subplot(1, 3, 1)
set_fig_properties([plt.gca()])
for runkey in hac.keys():
    plt.plot(hac[runkey]['t'], (hac[runkey]['kxi'] + hac[runkey]['kyi']) / 2, color='C7', alpha=0.5)
plt.plot(t, ki_ave, color='C3', linewidth=1)
plt.xlim([0, 1000])
plt.gca().set_xticks(range(0, 1001, 200))
plt.ylim([0, 5])
plt.xlabel('Correlation Time (ps)')
plt.ylabel(r'$\kappa^{in}$ (W/m/K)')
plt.title('(a)')

plt.subplot(1, 3, 2)
set_fig_properties([plt.gca()])
for runkey in hac.keys():
    plt.plot(hac[runkey]['t'], (hac[runkey]['kxo'] + hac[runkey]['kyo']) / 2, color='C7', alpha=0.5)
plt.plot(t, ko_ave, color='C0', linewidth=1)
plt.xlim([0, 1000])
plt.gca().set_xticks(range(0, 1001, 200))
plt.ylim([0, 5])
plt.xlabel('Correlation Time (ps)')
plt.ylabel(r'$\kappa^{out}$ (W/m/K)')
plt.title('(b)')

plt.subplot(1, 3, 3)
set_fig_properties([plt.gca()])
plt.plot(t, ko_ave, color='C0', linewidth=1)
plt.plot(t, ki_ave, color='C3', linewidth=1)
plt.plot(t, ki_ave + ko_ave, color='k', linewidth=1)
plt.xlim([0, 1000])
plt.gca().set_xticks(range(0, 1001, 200))
plt.ylim([0, 5])
plt.xlabel('Correlation Time (ps)')
plt.ylabel(r'$\kappa$ (W/m/K)')
plt.title('(c)')

plt.tight_layout()
plt.savefig('./EMD.png', bbox_inches='tight')

