# Usage: python latent_AB.py A.xyz B.xyz
import matplotlib.pyplot as plt
import numpy as np
import sys
from ase.io import read
from pynep.calculate import NEP
from sklearn.decomposition import PCA

plt.rcParams['font.size'] = 14
plt.rcParams["axes.unicode_minus"] = False
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'

def plt_latent(filename, p0, p1, proj, selected_proj):
    plt.figure(figsize=(5, 4), dpi=300)
    plt.scatter(proj[:,0], proj[:,1], c='blue', alpha=0.5, label='A')
    plt.scatter(selected_proj[:,0], selected_proj[:,1], c='orange', alpha=0.5, label='B')
    plt.xlabel('PCA dimension 0 - Var={:.2f}'.format(p0))
    plt.ylabel('PCA dimension 1 - Var={:.2f}'.format(p1))
    # plt.xlim(-1.10,1.10)
    # plt.ylim(-0.35,0.55)
    # plt.xticks(np.arange(-1.10, 1.15, 0.4))
    # plt.yticks(np.arange(-0.30, 0.55, 0.2))
    plt.legend()
    plt.tight_layout()
    plt.savefig('latent_{}.png'.format(filename), bbox_inches='tight')  

a = read(sys.argv[1], ':')
b = read(sys.argv[2], ':')
calc = NEP("./nep.txt")
print(calc)

des_a = np.array([np.mean(calc.get_property('descriptor', i), axis=0) for i in a])
des_b = np.array([np.mean(calc.get_property('descriptor', i), axis=0) for i in b])

reducer = PCA(n_components=2)
reducer.fit(des_a)
p0, p1 = reducer.explained_variance_ratio_[:2]
proj_a = reducer.transform(des_a)
proj_b = reducer.transform(des_b)
colors = [len(i) for i in a]
plt.scatter(proj_a[:,0], proj_a[:,1], c=colors, label='NEP train')
plt_latent('AB', p0, p1, proj_a, proj_b) 