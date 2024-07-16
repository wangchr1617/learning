import matplotlib.pyplot as plt
import numpy as np
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'

virial = np.loadtxt('./virial_train.out')
loss = np.loadtxt('./loss.out')
nep = virial[:, 0:6]
dft = virial[:, 6:12]

lim0, lim1 = int(dft.min())-1, int(dft.max())+1
distance = lim1 - lim0

plt.figure(figsize=(5,4))
plt.plot(dft, nep, 'o', markersize=3)
plt.plot(np.linspace(lim0, lim1), np.linspace(lim0, lim1), '-', c='#85618E')
plt.xlim(lim0, lim1)
plt.ylim(lim0, lim1)
plt.xlabel('DFT virial (eV/atom)')
plt.ylabel('NEP virial (eV/atom)')
string = 'train RMSE = {:.1f} meV/Å'.format(loss[-1,6]*1000)
plt.text(lim1-distance*0.05, lim0+distance*0.05, string, horizontalalignment='right', verticalalignment='baseline')
plt.legend(['xx', 'yy', 'zz', 'xy', 'yz', 'zx'])
plt.tight_layout()
plt.savefig('./virial.png')
