import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import norm
plt.rcParams['font.size'] = 14
plt.rcParams["axes.unicode_minus"] = False
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
import warnings
warnings.filterwarnings('ignore')

def force_err():
    r = np.loadtxt('./ref.txt').reshape((-1,1))
    d = np.loadtxt('./dft.txt').reshape((-1,1))
    n = np.loadtxt('./nep.txt').reshape((-1,1))
    n = np.concatenate((r, d, n), axis=1)
    df = pd.DataFrame(n, columns=["ref","dft","nep"])
    df["e_d"] = np.array(df["ref"]) - np.array(df["dft"])
    df["e_n"] = np.array(df["ref"]) - np.array(df["nep"])
    # df.to_csv('./data.csv')
    return df

def plot(df):
    e_d = np.array(df['e_d'])
    e_n = np.array(df['e_n'])
   
    plt.figure(figsize=(5, 4), dpi=300)
    # plt.axvline(x=0, color='k', linewidth=1, linestyle='-')    

    plt.hist(e_d, bins=30, color='k', alpha=0.3, density=True, label='DFT')
    xmin, xmax = plt.xlim()
    x_d = np.linspace(xmin, xmax, 100)
    mu_d, std_d = norm.fit(e_d)
    p_d = norm.pdf(x_d, mu_d, std_d)
    plt.plot(x_d, p_d, linewidth=1, linestyle='--', color='k', alpha=0.75)
    plt.fill_between(x_d, p_d, where=(p_d>=0), color='k', alpha=0.5)    
    
    plt.hist(e_n, bins=30, color='r', alpha=0.3, density=True, label='NEP')
    xmin, xmax = plt.xlim()
    x_n = np.linspace(xmin, xmax, 100)
    mu_n, std_n = norm.fit(e_n)
    p_n = norm.pdf(x_n, mu_n, std_n)
    plt.plot(x_n, p_n, linewidth=1, linestyle='--', color='r', alpha=0.75)
    plt.fill_between(x_n, p_n, where=(p_n>=0), color='r', alpha=0.5)    

    plt.xlabel('Force error (eV/A)')
    plt.ylabel('Density')
    plt.xlim(-0.25, 0.25)
    plt.xticks(np.arange(-0.2, 0.25, 0.1))
    plt.ylim(0, max(p_d))
    plt.legend()
    
    plt.tight_layout()
    plt.savefig('./Force_error.png', bbox_inches='tight')

df = force_err()
plot(df)

