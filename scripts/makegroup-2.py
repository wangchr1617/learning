# Usage: python *.py
def group_atoms(atoms, ranges):
    for atom in atoms:
        y = atom['z']
        for y_min, y_max, group_property in ranges:
            if y_min <= y < y_max:
                atom['property'] = group_property
                break

ranges = [(-1, 64, 0), (64, 127, 1)]

with open("model.xyz", 'r') as f:
    lines = f.readlines()

lattice_properties_line = lines[1].strip()
lattice_properties_line += ":group:I:1\n"

lines = lines[2:]
atoms = []
for line in lines:
    atom_data = line.split()
    atom = {
        'element': atom_data[0],
        'x': float(atom_data[1]),
        'y': float(atom_data[2]),
        'z': float(atom_data[3]),
        'property': None
    }
    atoms.append(atom)

group_atoms(atoms, ranges)

with open("model_new.xyz", 'w') as f:
    f.write(f"{len(atoms)}\n")
    f.write(lattice_properties_line)
    for atom in atoms:
        f.write(f"{atom['element']} {atom['x']} {atom['y']} {atom['z']} {atom['property']}\n")

print("New XYZ file has been generated.")
