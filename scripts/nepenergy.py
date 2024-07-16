import matplotlib.pyplot as plt
import numpy as np
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'

energy = np.loadtxt('./energy_train.out')
loss = np.loadtxt('./loss.out')
dft = energy[:, 1]
nep = energy[:, 0]

lim0, lim1 = int(dft.min())-1, int(dft.max())+1
distance = lim1 - lim0

plt.figure(figsize=(5,4))
plt.plot(dft, nep, 'o', c='#E57B71')
plt.plot(np.linspace(lim0, lim1), np.linspace(lim0, lim1), '-', c='#85618E')
plt.xlim(lim0, lim1)
plt.ylim(lim0, lim1)
plt.xlabel('DFT energy (eV/atom)')
plt.ylabel('NEP energy (eV/atom)')
string = 'train RMSE = {:.1f} meV/atom'.format(loss[-1,4]*1000)
plt.text(lim1-distance*0.05, lim0+distance*0.05, string, horizontalalignment='right', verticalalignment='baseline')
plt.tight_layout()
plt.savefig('./energy.png')
