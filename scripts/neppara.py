# Usage: python neppara.py
import matplotlib.pyplot as plt
import numpy as np
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'

para = np.loadtxt('./nep.restart')
x = np.arange(0, len(para))
y = para[:, 0]

plt.figure(figsize=(5,4))
plt.plot(x, y)
plt.xlim(min(x), max(x))
plt.tight_layout()
plt.savefig('./para.png')

