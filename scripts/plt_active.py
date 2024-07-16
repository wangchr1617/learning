import matplotlib.pyplot as plt
import numpy as np
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'

x,y = np.loadtxt('active.out', usecols=(0,1), unpack=True)
plt.xlabel('Idx')
plt.ylabel('Error')
plt.scatter(x, y)
plt.savefig('./active.png', bbox_inches='tight')
