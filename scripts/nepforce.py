import matplotlib.pyplot as plt
import numpy as np
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'

force = np.loadtxt('./force_train.out')
loss = np.loadtxt('./loss.out')
nep = force[:, 0:3]
dft = force[:, 3:6]

lim0, lim1 = int(dft.min())-10, int(dft.max())+10
distance = lim1 - lim0

plt.figure(figsize=(5,4))
plt.plot(dft, nep, 'o')
plt.plot(np.linspace(lim0, lim1), np.linspace(lim0, lim1), '-', c='#85618E')
plt.xlim(lim0, lim1)
plt.ylim(lim0, lim1)
plt.xlabel('DFT force (eV/Å)')
plt.ylabel('NEP force (eV/Å)')
string = 'train RMSE = {:.1f} meV/Å'.format(loss[-1,5]*1000)
plt.text(lim1-distance*0.05, lim0+distance*0.05, string, horizontalalignment='right', verticalalignment='baseline')
plt.legend(['x direction', 'y direction', 'z direction'])
plt.tight_layout()
plt.savefig('./force.png')
