import numpy as np

e = np.loadtxt('energy_train.out')
rmse_e = np.sqrt(np.mean((e[:,0] - e[:,1])**2))
print(rmse_e*1000, 'meV/atom')

f = np.loadtxt('force_train.out')
rmse_f = np.sqrt(np.mean((f[:,3:6] - f[:,0:3])**2))
print(rmse_f*1000, 'meV/A')

v = np.loadtxt('virial_train.out')
rmse_v = np.sqrt(np.mean((v[:,6:12] - v[:,0:6])**2))
print(rmse_v*1000, 'meV/atom')

