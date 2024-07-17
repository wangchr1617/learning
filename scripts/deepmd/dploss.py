import numpy as np
import matplotlib.pyplot as plt

data = np.genfromtxt("lcurve.out", names=True)
for name in data.dtype.names[1:-1]:
    plt.plot(data['step'], data[name], label=name)
plt.xlabel('Step')
plt.ylabel('Loss')
# plt.loglog()
plt.xscale('symlog')
plt.yscale('log')
plt.grid()
plt.legend()
plt.savefig('./loss.png', bbox_inches='tight')
