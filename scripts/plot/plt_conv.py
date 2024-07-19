# Usage: python *.py natom
# conv.txt: grep F= OSZICAR | awk '{print \$1,\$5}' > conv.txt

import matplotlib.pyplot as plt
import numpy as np
import sys
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'

data = np.loadtxt("./conv.txt", delimiter=' ', dtype=float)
natom = int(sys.argv[1])

plt.figure(figsize=(4,3))
plt.plot(data[:, 0], data[:, 1]/natom, c='r', alpha=0.8)
plt.scatter(data[:, 0], data[:, 1]/natom, c='b', s=3, alpha=0.8)
plt.xlim(data[0,0]-1, data[-1,0])
plt.xlabel('Steps')
plt.ylabel('Energy per atom')
plt.tight_layout()
plt.savefig('./conv.png')
