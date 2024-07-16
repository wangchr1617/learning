import numpy as np

def read_poscar_lattice(filename):
    lattice = np.zeros((3, 3))
    try:
        with open(filename, 'r') as file:
            file.readline()
            scale_factor = float(file.readline().strip())
            for i in range(3):
                line = file.readline().strip().split()
                lattice[i] = [float(x) * scale_factor for x in line[:3]]
    except IOError:
        print(f"Error: Could not read file {filename}")
    except ValueError:
        print("Error: File format is not correct")
    return lattice

unit = np.array(read_poscar_lattice("./POSCAR"))
print(unit)
prim = np.array(read_poscar_lattice("./PRIMCELL.vasp"))
print(prim)
unit_inv = np.linalg.inv(unit)
result = np.dot(unit_inv, prim)
print(result)
