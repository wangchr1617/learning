import numpy as np
import matplotlib.pyplot as plt
from ase.io import read, write
from pynep.calculate import NEP
from pynep.select import FarthestPointSample
from sklearn.decomposition import PCA

# Read data
a = read('./nepmd.xyz', ':')
atomsnum = [len(i) for i in a]
calc = NEP("./nep.txt")
print(calc)

# Calculate descriptors of atoms
descriptors = []
for i in range(len(a)):
    descriptors.append(calc.get_property('descriptor', a[i]))
all_descriptors = np.concatenate(descriptors, axis=0)
descriptor_mean = all_descriptors.mean(axis=0)
descriptor_std = all_descriptors.std(axis=0)

# Calculate descriptors of structures
all_descriptors_mean = np.array([np.mean(descriptors[i].reshape(1,-1), axis=0) for i in range(len(descriptors))])
normalized_descriptors_mean = np.array([np.mean(np.array([(d - descriptor_mean) / descriptor_std for d in descriptors[i]]).reshape(1,-1), axis=0) for i in range(len(descriptors))])

# Farthest point sample
def fps(des):
    sampler = FarthestPointSample(min_distance=0.05)
    selected_i = sampler.select(des, [], max_select=100)
    unselected_i = []
    for i in range(len(a)):
        if i not in selected_i:
            unselected_i.append(i)
    #write('./selected.xyz', [a[i] for  i in selected_i])
    #write('./unselected.xyz', [a[i] for i in unselected_i])
    return selected_i, unselected_i

# PCA
fig, axes = plt.subplots(ncols=2, figsize=(7, 3), dpi=140)
for k, ax in enumerate(axes):
    pca = PCA(n_components=2)
    if k == 0:
        pc = pca.fit_transform(all_descriptors_mean)
        selected_i, _ = fps(all_descriptors_mean)
        selected_pc = pca.transform(np.array([all_descriptors_mean[i] for i in selected_i]))
        title = 'unnormalized'
    else:
        pc = pca.fit_transform(normalized_descriptors_mean)
        selected_i, _ = fps(normalized_descriptors_mean)
        selected_pc = pca.transform(np.array([normalized_descriptors_mean[i] for i in selected_i]))
        title = 'normalized'
    ax.scatter(pc[:, 0], pc[:, 1], s=0.8, c='b', alpha=0.5, label='raw data')
    ax.scatter(selected_pc[:, 0], selected_pc[:, 1], s=0.8, c='r', alpha=0.5, label='selected data')
    ax.set_xlabel(f'PCA dimension 0 - Var={pca.explained_variance_ratio_[0]:.2f}')
    ax.set_ylabel(f'PCA dimension 1 - Var={pca.explained_variance_ratio_[1]:.2f}')
    ax.set_title(title)
ax.legend(frameon=False)
fig.align_ylabels()
#plt.colorbar()
plt.tight_layout()
plt.savefig('fps.png')        
