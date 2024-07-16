# Usage: python bondlengths_distribution.py 100.xyz 100
from ase.data import covalent_radii
from ase.io import read
from ase.neighborlist import NeighborList

import matplotlib.pyplot as plt
import sys

def get_bondpairs(atoms, radius=1.1):
    cutoffs = radius * covalent_radii[atoms.numbers]
    nl = NeighborList(cutoffs=cutoffs, self_interaction=False)
    nl.update(atoms)
    positions = atoms.get_positions()
    bondpairs = []
    for a in range(len(atoms)):
        indices, _ = nl.get_neighbors(a)
        for b in indices:
            dis = atoms.get_distance(a, b, mic=True)
            bondpairs.append((a, b, dis))
    return bondpairs

traj = read(sys.argv[1], index=':')
idx = sys.argv[2]
radius = 1.1
bondlengths = []
for atoms in traj:
    bondpairs = get_bondpairs(atoms, radius)
    bondlengths = bondlengths + [pair[2] for pair in bondpairs]

fig, ax = plt.subplots()
ax.hist(bondlengths, bins=50, density=True, alpha=0.5, color='r', edgecolor='black')
ax.set_xlabel('Bond Length')
ax.set_ylabel('Frequency')
plt.savefig('./bondlengths_{}.png'.format(idx))
