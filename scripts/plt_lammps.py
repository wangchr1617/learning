import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from matplotlib.ticker import MultipleLocator, FormatStrFormatter
plt.rcParams['font.size'] = 10
plt.rcParams["axes.unicode_minus"] = False
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'

colors = sns.color_palette("Paired")[4:10]

file_name = './log.lammps'
natoms = 216
data = []
cols = ["Step", "PotEng", "KinEng", "TotEng", "Temp", "Press", "Volume"]

with open(file_name, 'r') as f:
    flag = 0
    for line in f:
        if "Step         PotEng         KinEng         TotEng          Temp          Press          Volume" in line:
            flag = 1
            continue
        if "Loop time" not in line and flag == 1:
            data.append(line.split())  
        else:
            flag = 0
df = pd.DataFrame(data, columns=cols, dtype='float')
df["Step"] = np.array(df["Step"]) * 0.5 * 0.001
df["PotEng"] = np.array(df["PotEng"]) / natoms
df["KinEng"] = np.array(df["KinEng"]) / natoms
df["TotEng"] = np.array(df["TotEng"]) / natoms

df.drop(df.index[[61, 122, 263, 324, 505]], inplace=True)
df.reset_index(drop=True, inplace=True)
df = df[:801]

plt.figure(dpi=600)
plt.plot(df["Step"][:62], df["PotEng"][:62], c=colors[0], label="2000 K")
plt.plot(df["Step"][61:121], df["PotEng"][61:121], c=colors[1], label="1000 K")
plt.plot(df["Step"][121:261], df["PotEng"][121:261], c=colors[2], label="Quenching")
plt.plot(df["Step"][261:321], df["PotEng"][261:321], c=colors[3], label="300 K")
plt.plot(df["Step"][321:501], df["PotEng"][321:501], c=colors[4], label="Annealing")
plt.plot(df["Step"][501:], df["PotEng"][501:], c=colors[5], label="600K")

plt.xlabel('Time (ps)')
plt.ylabel('Potential Energy (eV/atom)')
plt.xlim(0, 400)
plt.gca().yaxis.set_major_locator(MultipleLocator(0.1))
plt.gca().yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
plt.legend()
plt.savefig('./PotEne-Time.png', bbox_inches='tight')
