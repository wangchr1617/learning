# Usage: python convergence_plot.py ./opt/ (POSCAR_ini)

from ase.io import read
from ase.io.vasp import read_vasp, read_vasp_xdatcar
import numpy as np
import matplotlib.pyplot as plt
import sys
import warnings
warnings.filterwarnings('ignore')
plt.rcParams['font.size'] = 14
# plt.rcParams['font.weight'] = 'bold'
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'

def WriteForce(outdata, atoms):
    fmax_list = []
    f = []
    for x, line in enumerate(outdata):
        if "TOTAL-FORCE" in line:
            f = []
            for y in range(x+2, len(outdata)):
                if outdata[y][1] == "-":
                    break
                f.append(outdata[y].split()[3:6])
            try:
                c = atoms._get_constraints()
                indices_fixed = c[0].index
                for i in indices_fixed:
                    f[i] = [0,0,0]
            except:
                pass
            fmax = 0
            for i in f:
                fval = (float(i[0])**2 + float(i[1])**2 + float(i[2])**2)**(1./2.)
                if fval > fmax:
                    fmax = fval
            fmax_list.append(fmax)
    return fmax_list

def Cartesian2Fractional(structure):
    cell = np.matrix(structure.cell)
    X = (cell.T)**-1
    cart_coords = np.matrix(structure.get_positions())
    frac_coords = (X*cart_coords.T).T
    return np.array(frac_coords)

path = sys.argv[1]

# Energy convergence
ion_steps = []
energies = []
with open("{}OSZICAR".format(path), "r") as f:
    for line in f.readlines():
        if "E0" in line:
            parts = line.split()
            ion_steps.append(int(parts[0]))
            energies.append(float(parts[4]))
plt.figure(figsize=(10, 6))
plt.grid(True)
plt.plot(ion_steps, energies, '-o', label='Energy')
plt.title('Energy of each ion steps')
plt.xlabel('Ion steps')
plt.ylabel('Energy(eV)')
plt.xlim(0, max(ion_steps)+1)
plt.legend()
plt.savefig('conv_ene.png')

# Force convergence
out = open("{}OUTCAR".format(path)).readlines()
atoms = read("{}CONTCAR".format(path))
forces = WriteForce(out, atoms)
plt.figure(figsize=(10, 6))
plt.grid(True)
plt.plot(range(1, len(forces)+1), forces, "-o", label="Force")
plt.title('Max Force vs. Optimization Step')
plt.xlabel('Optimization Step')
plt.ylabel('Max Force (eV/Å)')
plt.xlim(0, len(forces)+1) 
plt.legend()
plt.savefig('conv_force.png')

# Distance from initial structure
if len(sys.argv) == 2:
    initial_structure = read_vasp("{}POSCAR".format(path))
elif len(sys.argv) == 3:
    print("Warning! POSCAR is not the initial structure you difined.")
    initial_structure = read_vasp(sys.argv[2])
structures = read_vasp_xdatcar("{}XDATCAR".format(path), index=slice(None))
distances = []
for structure in structures:
    rmsd = np.sqrt(((Cartesian2Fractional(structure) - Cartesian2Fractional(initial_structure)) ** 2).sum(axis=1).mean())
    distances.append(rmsd)
plt.figure(figsize=(10, 6))
plt.grid(True)
plt.plot(range(1, len(distances)+1), distances, "-o", label='Structural change')
plt.title('Distance from initial structure')
plt.xlabel('Optimization step')
plt.ylabel('RMSD from initial structure')
plt.xlim(0, len(distances)+1) 
plt.legend()
plt.savefig('conv_structure.png')

print("Over!!!")
