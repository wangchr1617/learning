
from gpyumd.load import load_shc, load_kappa
from gpyumd.math import running_ave
from gpyumd.calc import calc_spectral_kappa

import numpy as np
import matplotlib.pyplot as plt

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

kappa = load_kappa()
kappa.keys()

t = np.arange(1, kappa['kxi'].shape[0] + 1) * 0.001  # ns
kappa['kyi_ra'] = running_ave(kappa['kyi'], t)
kappa['kyo_ra'] = running_ave(kappa['kyo'], t)
kappa['kxi_ra'] = running_ave(kappa['kxi'], t)
kappa['kxo_ra'] = running_ave(kappa['kxo'], t)
kappa['kz_ra'] = running_ave(kappa['kz'], t)

plt.figure(figsize=(12, 10))
plt.subplot(2, 2, 1)
set_fig_properties([plt.gca()])
plt.plot(t, kappa['kyi'], color='C7', alpha=0.5)
plt.plot(t, kappa['kyi_ra'], linewidth=2)
plt.xlim([0, 10])
plt.gca().set_xticks(range(0, 11, 2))
plt.ylim([-10, 10])
# plt.gca().set_yticks(range(-2000, 4001, 1000))
plt.xlabel('time (ns)')
plt.ylabel(r'$\kappa_{in}$ W/m/K')
plt.title('(a)')

plt.subplot(2, 2, 2)
set_fig_properties([plt.gca()])
plt.plot(t, kappa['kyo'], color='C7', alpha=0.5)
plt.plot(t, kappa['kyo_ra'], linewidth=2, color='C3')
plt.xlim([0, 10])
plt.gca().set_xticks(range(0, 11, 2))
plt.ylim([-10, 10])
# plt.gca().set_yticks(range(0, 4001, 1000))
plt.xlabel('time (ns)')
plt.ylabel(r'$\kappa_{out}$ (W/m/K)')
plt.title('(b)')

plt.subplot(2, 2, 3)
set_fig_properties([plt.gca()])
plt.plot(t, kappa['kyi_ra'], linewidth=2)
plt.plot(t, kappa['kyo_ra'], linewidth=2, color='C3')
plt.plot(t, kappa['kyi_ra'] + kappa['kyo_ra'], linewidth=2, color='k')
plt.xlim([0, 10])
plt.gca().set_xticks(range(0, 11, 2))
plt.ylim([-10, 10])
# plt.gca().set_yticks(range(0, 4001, 1000))
plt.xlabel('time (ns)')
plt.ylabel(r'$\kappa$ (W/m/K)')
plt.legend(['in', 'out', 'total'])
plt.title('(c)')

plt.subplot(2, 2, 4)
set_fig_properties([plt.gca()])
plt.plot(t, kappa['kyi_ra'] + kappa['kyo_ra'], color='k', linewidth=2)
plt.plot(t, kappa['kxi_ra'] + kappa['kxo_ra'], color='C0', linewidth=2)
plt.plot(t, kappa['kz_ra'], color='C3', linewidth=2)
plt.xlim([0, 10])
plt.gca().set_xticks(range(0, 11, 2))
plt.ylim([-10, 10])
# plt.gca().set_yticks(range(-2000, 4001, 1000))
plt.xlabel('time (ns)')
plt.ylabel(r'$\kappa$ (W/m/K)')
plt.legend(['yy', 'xy', 'zy'])
plt.title('(d)')

plt.tight_layout()
plt.savefig('./HNEMD.png', bbox_inches='tight')

