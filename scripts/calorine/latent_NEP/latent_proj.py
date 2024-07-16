# Usage: python latent_proj.py A.xyz B.xyz
import numpy as np
import pandas as pd
import sys
from ase.io import read
from pynep.calculate import NEP
from sklearn.decomposition import PCA

a = read(sys.argv[1], ':')
b = read(sys.argv[2], ':')
calc = NEP("./nep.txt")
print(calc)

des_a = np.array([np.mean(calc.get_property('descriptor', i), axis=0) for i in a])
des_b = np.array([np.mean(calc.get_property('descriptor', i), axis=0) for i in b])

reducer = PCA(n_components=2)
reducer.fit(des_a)
p0, p1 = reducer.explained_variance_ratio_[:2]
print(p0, p1)
proj_b = reducer.transform(des_b)
df_b = pd.DataFrame(proj_b[:,0:2])
df_b.to_csv("proj_b.csv")
