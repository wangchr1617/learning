import matplotlib.pyplot as plt
import numpy as np
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'

cols = ['Generation', 'Total', 'L1-regularization', 'L2-regularization', 'Energy-train', 'Force-train',
        'Virial-train', 'Energy-test', 'Force-test', 'Virial-test']
loss = np.loadtxt('loss.out')

plt.figure(figsize=(5,4))
plt.plot(loss[:, 0], loss[:, 1:7])
plt.loglog()
plt.xlim(100, loss[-1, 0])
plt.xlabel('Generation')
plt.ylabel('Loss')
plt.legend(cols[1:7])
plt.tight_layout()
plt.savefig('./loss.png')
