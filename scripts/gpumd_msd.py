import numpy as np
import matplotlib.pyplot as plt

plt.rcParams['font.size'] = 10
plt.rcParams["axes.unicode_minus"] = False
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'

file_path = 'msd.out'

data = np.loadtxt(file_path)

correlation_time = data[:, 0]
msd_x = data[:, 1]
msd_y = data[:, 2]
msd_z = data[:, 3]

plt.figure(figsize=(10, 6))
plt.plot(correlation_time, msd_x, label='MSD in x direction')
plt.plot(correlation_time, msd_y, label='MSD in y direction')
plt.plot(correlation_time, msd_z, label='MSD in z direction')

plt.xlabel('Correlation Time (ps)')
plt.ylabel('MSD (Å²)')
plt.title('MSD vs Correlation Time')
plt.legend()
plt.xlim(0, ) 
plt.ylim(0, ) 
plt.savefig("msd.png", bbox_inches='tight')

