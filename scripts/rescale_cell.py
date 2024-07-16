# Usage: python rescale_cell.py
from ase.io import read, write

poscar_file = 'POSCAR'
original_structure = read(poscar_file)

scaling_factors = [0.9, 0.92, 0.94, 0.96, 0.98, 1.0]

for factor in scaling_factors:
    modified_structure = original_structure.copy()
    cell = modified_structure.cell
    new_cell = cell.copy()
    new_cell[2, 2] *= factor
    modified_structure.set_cell(new_cell, scale_atoms=True)
    output_filename = f"{factor:.2f}.vasp"
    write(output_filename, modified_structure, format='vasp', direct=True, vasp5=True)

print("VASP files generated successfully.")

