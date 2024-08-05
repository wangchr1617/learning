import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

plt.rcParams['font.size'] = 10
plt.rcParams["axes.unicode_minus"] = False
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'

palette = [[0.96779756, 0.44127456, 0.53581032],
           [0.80879541, 0.56347001, 0.19502643],
           [0.59208915, 0.6418467 , 0.19350691],
           [0.19783576, 0.6955517 , 0.3995301 ],
           [0.21044754, 0.67731051, 0.64339412],
           [0.22335772, 0.65657923, 0.81713555],
           [0.64230443, 0.54976801, 0.95826514],
           [0.96038885, 0.38143179, 0.86831177]]

def volume(a, b, c, ncell):
    arr = np.array([a, b, c]).T
    det = np.linalg.det(arr) / ncell
    return det

def length(vec, ncell):
    return np.linalg.norm(vec) / ncell

def angle(vec1, vec2):
    cos_theta = np.dot(vec1, vec2) / (np.linalg.norm(vec1) * np.linalg.norm(vec2))
    angle_radians = np.arccos(np.clip(cos_theta, -1.0, 1.0))
    return np.abs(90 - np.degrees(angle_radians))

def plot_md(filename):
    data = np.loadtxt(filename)
    cols = ["T", "K", "U", "Px", "Py", "Pz", "Pyz", "Pxz", "Pxy", 
            "ax", "ay", "az", "bx", "by", "bz", "cx", "cy", "cz"]
    thermo = pd.DataFrame(data, columns=cols)
    
    ncell = 12
    natom = (ncell**3) * 8
    temp = np.arange(100, 800, 1)[::-1]
    interval = int(len(thermo) / len(temp))
    thermo = thermo[::interval].reset_index(drop=True)
    
    v_ave = np.array([volume(thermo.loc[i, ['ax', 'ay', 'az']], 
                             thermo.loc[i, ['bx', 'by', 'bz']], 
                             thermo.loc[i, ['cx', 'cy', 'cz']], ncell**3) 
                      for i in range(len(thermo))])
    
    a_ave = np.array([length(thermo.loc[i, ['ax', 'ay', 'az']], ncell) for i in range(len(thermo))])
    b_ave = np.array([length(thermo.loc[i, ['bx', 'by', 'bz']], ncell) for i in range(len(thermo))])
    c_ave = np.array([length(thermo.loc[i, ['cx', 'cy', 'cz']], ncell) for i in range(len(thermo))])
    
    alpha = np.array([angle(thermo.loc[i, ['bx', 'by', 'bz']], thermo.loc[i, ['cx', 'cy', 'cz']]) 
                      for i in range(len(thermo))])
    beta = np.array([angle(thermo.loc[i, ['ax', 'ay', 'az']], thermo.loc[i, ['cx', 'cy', 'cz']]) 
                     for i in range(len(thermo))])
    gamma = np.array([angle(thermo.loc[i, ['ax', 'ay', 'az']], thermo.loc[i, ['bx', 'by', 'bz']]) 
                      for i in range(len(thermo))])
    
    p_ave = (thermo[['Px', 'Py', 'Pz', 'Pyz', 'Pxz', 'Pxy']].sum(axis=1)) / 6.0
    px = np.array(thermo['Px'])
    py = np.array(thermo['Py'])
    pz = np.array(thermo['Pz'])
    pyz = np.array(thermo['Pyz'])
    pxz = np.array(thermo['Pxz'])
    pxy = np.array(thermo['Pxy'])
    
    tot_ene = (thermo['K'] + thermo['U']) / natom
    
    plt.figure(figsize=(10, 8), dpi=300)
    
    plt.subplot(2, 2, 1)
    plt.plot(temp, v_ave, c='k')
    plt.axvline(275, c='k', linestyle='--', linewidth=1, alpha=0.5)
    plt.axvline(450, c='k', linestyle='--', linewidth=1, alpha=0.5)
    plt.xlim(100, 800)
    # plt.ylim(224, 229)
    plt.xlabel('Temperature (K)')
    plt.ylabel(r'Normalized Lattice Volume (Å³/f.u.)')
    plt.xticks(np.arange(100, 850, 100))
    plt.title('(a)')
    
    plt.subplot(2, 2, 2)
    plt.plot(temp, a_ave, c=palette[1], label='a')
    plt.plot(temp, b_ave, c=palette[2], label='b')
    plt.plot(temp, c_ave, c=palette[3], label='c')
    plt.axvline(275, c='k', linestyle='--', linewidth=1, alpha=0.5)
    plt.axvline(450, c='k', linestyle='--', linewidth=1, alpha=0.5)
    plt.xlim(100, 800)
    # plt.ylim(6.05, 6.13)
    plt.xlabel('Temperature (K)')
    plt.ylabel(r'Normalized Lattice Length (Å)')
    plt.xticks(np.arange(100, 850, 100))
    plt.legend(loc="lower right")
    plt.title('(b)')
    
    plt.subplot(2, 2, 3)
    plt.plot(temp, alpha, c=palette[1], label=r'$\alpha$')
    plt.plot(temp, beta,  c=palette[2], label=r'$\beta$')
    plt.plot(temp, gamma, c=palette[3], label=r'$\gamma$')
    plt.axvline(275, c='k', linestyle='--', linewidth=1, alpha=0.5)
    plt.axvline(450, c='k', linestyle='--', linewidth=1, alpha=0.5)
    plt.xlim(100, 800)
    plt.ylim(-0.1, )
    plt.xlabel('Temperature (K)')
    plt.ylabel(r'90°- Lattice Angle (°)')
    plt.xticks(np.arange(100, 850, 100))
    plt.legend(loc="upper right")
    plt.title('(c)')
    
    plt.subplot(2, 2, 4)
    plt.plot(temp, px, c=palette[0], alpha=0.5, label='Px')
    plt.plot(temp, py, c=palette[1], alpha=0.5, label='Py')
    plt.plot(temp, pz, c=palette[2], alpha=0.5, label='Pz')
    plt.plot(temp, pyz, c=palette[3], alpha=0.5, label='Pyz')
    plt.plot(temp, pxz, c=palette[4], alpha=0.5, label='Pxz')
    plt.plot(temp, pxy, c=palette[5], alpha=0.5, label='Pxy')
    plt.plot(temp, p_ave, c='k', label='P_ave')
    plt.xlim(100, 800)
    plt.ylim(-0.5, 0.5)
    plt.xlabel('Temperature (K)')
    plt.ylabel('Pressure (GPa)')
    plt.xticks(np.arange(100, 850, 100))
    plt.legend(loc="upper right")
    plt.title('(d)')
    
    plt.tight_layout()
    plt.savefig('./expansion.png', bbox_inches='tight')

plot_md('./thermo.out')
