# Usage: python calculate_para_diffs.py NEP-dataset.xyz 56

import numpy as np
import argparse

def read_xyz(filename, n_atoms):
    with open(filename, 'r') as file:
        lines = file.readlines()

    frames = []
    energies = []
    virials = []
    current_frame = []
    current_energy = None
    current_virial = None

    for i, line in enumerate(lines):
        if line.strip().isdigit():
            if current_frame:
                frames.append(current_frame)
                energies.append(current_energy)
                virials.append(current_virial)
                current_frame = []
                current_energy = None
                current_virial = None
        elif "Lattice=" in line:
            parts = line.split()
            try:
                # Extract energy, ignore free_energy
                current_energy = float([p.split('=')[1] for p in parts if p.startswith("energy=")][0])
                # Extract virial and handle double quotes
                virial_str = [p.split('=')[1] for p in parts if p.startswith("virial=")][0]
                virial_str = virial_str.strip('"')
                current_virial = [float(v) for v in virial_str.split()]
            except (IndexError, ValueError) as e:
                print(f"Error parsing line {i}: {line}")
                print(f"Parts: {parts}")
                raise e
        elif len(line.split()) == 7:  # This line contains atom data
            current_frame.append(list(map(float, line.split()[1:])))

    if current_frame:
        frames.append(current_frame)
        energies.append(current_energy)
        virials.append(current_virial)

    print(f"Total frames read: {len(frames)}")
    return np.array(frames), np.array(energies), np.array(virials)

def calculate_rmse(frames, reference_frame):
    reference_forces = frames[reference_frame]
    n_frames = len(frames)

    forces_rmse = []
    for i in range(n_frames):
        if i == reference_frame:
            forces_rmse.append((i, 0.0))
            continue
        diff = frames[i] - reference_forces
        mse = np.mean(np.square(diff))
        rmse = np.sqrt(mse)
        forces_rmse.append((i, rmse))

    return forces_rmse

def calculate_energy_differences(reference_energy_per_atom, energies, n_atoms):
    differences = np.abs((energies / n_atoms) - reference_energy_per_atom)
    return list(enumerate(differences))

def calculate_virial_differences(reference_virial_per_atom, virials, n_atoms):
    differences = [np.mean(np.abs((np.array(v) / n_atoms) - reference_virial_per_atom)) for v in virials]
    return list(enumerate(differences))

def write_to_file(filename, data):
    with open(filename, 'w') as f:
        f.write("Index Value\n")
        for idx, value in data:
            f.write(f"{idx} {value:.6f}\n")

def main(filename, n_atoms):
    frames, energies, virials = read_xyz(filename, n_atoms)
    reference_frame = len(frames) - 1

    forces_rmse = calculate_rmse(frames, reference_frame)

    reference_energy_per_atom = energies[reference_frame] / n_atoms
    reference_virial_per_atom = np.array(virials[reference_frame]) / n_atoms

    energy_diffs = calculate_energy_differences(reference_energy_per_atom, energies, n_atoms)
    virial_diffs = calculate_virial_differences(reference_virial_per_atom, virials, n_atoms)

    write_to_file('forces.txt', forces_rmse)
    write_to_file('energy.txt', energy_diffs)
    write_to_file('virial.txt', virial_diffs)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Calculate RMSE, Energy/Atom Difference, and Virial/Atom Difference for XYZ frames.')
    parser.add_argument('filename', type=str, help='The XYZ file name')
    parser.add_argument('n_atoms', type=int, help='Number of atoms per frame')

    args = parser.parse_args()

    main(args.filename, args.n_atoms)

