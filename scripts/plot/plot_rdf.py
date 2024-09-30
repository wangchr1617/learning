# Usage: python plt_rdf.py
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'

data = np.loadtxt('./PCF.dat', skiprows=1, dtype=float)
cols = ["Distance", "RDF"]
df = pd.DataFrame(data, columns=cols)

plt.figure(figsize = (5,3))
plt.plot(df[cols[0]], df[cols[1]], label="RDF")
plt.xlim(0, max(df[cols[0]]))
plt.ylim(0, max(df[cols[1]])+1)
plt.xlabel("Distance / Angstrom")
plt.ylabel("RDF")
plt.title("Radial distribution function")
plt.legend()
plt.savefig("./rdf.png", bbox_inches='tight')
