# Usage: python latent_fps.py [s|a|cf|ct]
import math
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import numpy as np
import re
import sys
from ase.io import read
from pynep.select import FarthestPointSample
from pynep.calculate import NEP
from sklearn.decomposition import PCA

plt.rcParams['font.size'] = 14
plt.rcParams["axes.unicode_minus"] = False
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'

def read_xyz(filename):
    with open(filename) as fileobj:
        lines = fileobj.readlines()
    data = []
    idx = 0
    while len(lines) > 0:
        images = {}
        images["idx"] = idx
        natoms = int(lines.pop(0))
        images["natoms"] = natoms
        comment = lines.pop(0)
        images["comment"] = comment        
        config_type_match = re.search(r'config_type=([^\s]+)', comment)
        if config_type_match:
            config_type_value = config_type_match.group(1)
            images["config_type"] = config_type_value
        energy_match = re.search(r'energy=(-?\d+\.\d+)', comment)
        if energy_match:
            energy_value = eval(energy_match.group(1)) / natoms
            images["energy"] = energy_value        
        atoms = []
        for _ in range(images["natoms"]):
            line = lines.pop(0)
            atoms.append(line)
        images["atoms"] = atoms
        data.append(images)
        idx += 1
    return data

def write_xyz(filename, data):
    with open(filename, 'w') as fileobj:
        for images in data:
            natoms = images["natoms"]
            comment = images["comment"].rstrip()
            if '\n' in comment:
                raise ValueError('Comment line should not have line breaks.')
            fileobj.write('%d\n%s\n' % (natoms, comment))
            for atoms in images["atoms"]:
                line = atoms.rstrip()
                fileobj.write('%s\n' % (line))

def plt_latent(filename, p0, p1, proj, selected_proj):
    plt.figure(figsize=(5, 4), dpi=300)
    plt.scatter(proj[:,0], proj[:,1], c='blue', alpha=0.5, label='all')
    plt.scatter(selected_proj[:,0], selected_proj[:,1], c='orange', alpha=0.5, label='selected')
    plt.xlabel('PCA dimension 0 - Var={:.2f}'.format(p0))
    plt.ylabel('PCA dimension 1 - Var={:.2f}'.format(p1))
    # plt.xlim(-1.10,1.10)
    # plt.ylim(-0.35,0.55)
    # plt.xticks(np.arange(-1.10, 1.15, 0.4))
    # plt.yticks(np.arange(-0.30, 0.55, 0.2))
    plt.legend()
    plt.tight_layout()
    plt.savefig('latent_{}.png'.format(filename), bbox_inches='tight')  

def latent_structure(filename, calc, min_select):
    frames = read(filename, ':')
    des = np.array([np.mean(calc.get_property('descriptor', atoms), axis=0) for atoms in frames])
    sampler = FarthestPointSample()
    selected = sampler.select(des, [], min_select=math.floor(len(des) * min_select))
    unselected = [i for i in range(len(frames)) if i not in selected]
    reducer = PCA(n_components=2)
    reducer.fit(des)
    p0, p1 = reducer.explained_variance_ratio_[:2]
    proj = reducer.transform(des)
    selected_proj = reducer.transform(np.array([des[i] for i in selected]))
    plt_latent('structure', p0, p1, proj, selected_proj)  
    data = read_xyz(filename)
    write_xyz('./selected.xyz', [data[i] for i in selected])
    write_xyz('./unselected.xyz', [data[i] for i in unselected])

def latent_atomic(filename, calc, min_select):
    frames = read(filename, ':')
    des = np.concatenate([calc.get_property('latent', atoms) for atoms in frames])
    structure = np.concatenate([[i] * len(atoms) for i, atoms in enumerate(frames)])
    sampler = FarthestPointSample()
    selected = set([structure[i] for i in sampler.select(des, [], min_select=math.floor(len(des) * min_select))])
    unselected = [i for i in range(len(frames)) if i not in selected]
    reducer = PCA(n_components=2)
    reducer.fit(des)
    p0, p1 = reducer.explained_variance_ratio_[:2]
    proj = reducer.transform(des)
    indices = [i for i, s in enumerate(structure) if s in selected]
    selected_proj = reducer.transform(np.array([des[i] for i in indices]))
    plt_latent('atomic', p0, p1, proj, selected_proj)
    data = read_xyz(filename)
    write_xyz('./selected.xyz', [data[i] for i in selected])
    write_xyz('./unselected.xyz', [data[i] for i in unselected])
    
def latent_configtype(filename, calc):
    data = read_xyz(filename)
    config_types = set(images["config_type"] for images in data)
    dic = {config_type: {"energy": [], "idx": []} for config_type in config_types}
    for images in data:
        config_type = images["config_type"]
        idx = images["idx"]
        dic[config_type]["idx"].append(idx)
    frames = read(filename, ':')
    des = np.array([np.mean(calc.get_property('descriptor', atoms), axis=0) for atoms in frames])
    reducer = PCA(n_components=2)
    reducer.fit(des)
    p0, p1 = reducer.explained_variance_ratio_[:2]    
    plt.figure(figsize=(5, 4), dpi=300)
    for config_type in dic:
        selected = dic[config_type]["idx"]
        selected_proj = reducer.transform(np.array([des[i] for i in selected]))
        plt.scatter(selected_proj[:,0], selected_proj[:,1], alpha=0.5, label=config_type.lower())
    plt.xlabel('PCA dimension 0 - Var={:.2f}'.format(p0))
    plt.ylabel('PCA dimension 1 - Var={:.2f}'.format(p1))
    plt.xlim(-0.5, 0.5)
    plt.ylim(-0.25, 0.25)
    plt.xticks(np.arange(-0.5, 0.51, 0.2))
    plt.yticks(np.arange(-0.2, 0.25, 0.1))
    plt.legend(fontsize=12)
    plt.tight_layout()
    plt.savefig('latent_configtype.png', bbox_inches='tight')
    
def latent_contour(filename, calc):    
    frames = read(filename, ':')
    des = np.array([np.mean(calc.get_property('descriptor', atoms), axis=0) for atoms in frames])
    reducer = PCA(n_components=2)
    reducer.fit(des)
    p0, p1 = reducer.explained_variance_ratio_[:2]
    proj = reducer.transform(des)
    x = proj[:,0]
    y = proj[:,1]
    data = read_xyz(filename)
    z = np.array([images["energy"] for images in data])
    ngrid = 200
    xi = np.linspace(-1, 1, ngrid)
    yi = np.linspace(-1, 1, ngrid)
    triang = tri.Triangulation(x, y)
    interpolator = tri.LinearTriInterpolator(triang, z)
    # interpolator = tri.CubicTriInterpolator(triang, z)
    Xi, Yi = np.meshgrid(xi, yi)
    zi = interpolator(Xi, Yi)
    plt.figure(figsize=(5, 4), dpi=300)
    cntrf=plt.contourf(xi, yi, zi, levels=15, alpha=0.75, cmap=plt.cm.hot)
    plt.colorbar(cntrf)
    plt.plot(x, y, 'ko', ms=1)
    plt.xlabel('PCA dimension 0 - Var={:.2f}'.format(p0))
    plt.ylabel('PCA dimension 1 - Var={:.2f}'.format(p1))
    plt.xlim(-0.3, 0.1)
    plt.ylim(-0.2, 0.2)
    plt.xticks(np.arange(-0.3, 0.11, 0.1))
    plt.yticks(np.arange(-0.2, 0.21, 0.1))
    plt.tight_layout()
    plt.savefig('latent_contour.png', bbox_inches='tight')

if len(sys.argv) < 2 or sys.argv[1] not in ["s", "a", "cf", "ct"]:
    print("Usage: python phonbycalorine.py [s|a|cf|ct]")
    sys.exit(1)

latent_type = sys.argv[1]
filename = './train.xyz'
calc = NEP("./nep.txt")
min_select = 0.8
if latent_type == "s":
    latent_structure(filename, calc, min_select)
elif latent_type == "a":
    latent_atomic(filename, calc, min_select)
elif latent_type == "cf":
    latent_configtype(filename, calc)
elif latent_type == "ct":
    latent_contour(filename, calc)
