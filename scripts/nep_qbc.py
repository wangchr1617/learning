from ase.io import read, write
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
plt.rcParams['font.size'] = 14
plt.rcParams["axes.unicode_minus"] = False
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
import warnings
warnings.filterwarnings('ignore')

def qbc():
    n1 = np.loadtxt('./nep_energy-1.txt', usecols=(1)).reshape((-1, 1))
    n4 = np.loadtxt('./dft_energy.txt', usecols=(0)).reshape((-1, 1))
    n5 = np.loadtxt('./atoms_num.txt', usecols=(0)).reshape((-1, 1))
    n = np.concatenate((n1, n4, n5), axis=1)
    df = pd.DataFrame(n, columns=["nep1", "dft", "num"])
    with open('./config_type.txt', 'r') as file:
        config_type = file.read().splitlines()
    df["config_type"] = config_type[:len(df)]
    df["dft"] = np.array(df["dft"]) / np.array(df["num"])
    df["std"] = df[["nep1", "dft"]].std(axis=1, numeric_only=True)
    df = df.sort_values(by=["std"], ascending=[False])
    df.to_csv('./std.csv', index=False)
    return df

def read_df():
  df = pd.read_csv('std.csv', index_col=0)
  return df

def select(df, begin, end=None):
  df_index_list = df.index.tolist()[begin:end]
  a = read('./model.xyz', ':')
  write('qbc.xyz', [a[i] for  i in df_index_list])

def plot(df, begin=None, end=None):
  plt.figure(figsize=(5, 4), dpi=300)
  x = df['dft'][begin:end]
  y = df['nep1'][begin:end]
  xmax = max([max(x), max(y)]) * 0.995
  xmin = min([min(x), min(y)]) * 1.005
  plt.plot([xmin,xmax], [xmin,xmax], linewidth=1)
  plt.scatter(x, y, marker='o', s=1.5)
  plt.xlabel('DFT energy (eV/atom)')
  plt.ylabel('NEP energy (eV/atom)')
  plt.xlim(xmin, xmax)
  plt.ylim(xmin, xmax)
  plt.savefig('./energy.png', bbox_inches='tight')

df = qbc()
# df = read_df()
# select(df, begin=50)
plot(df)



