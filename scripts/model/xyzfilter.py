def read_xyz(file_path):
    frames = []
    with open(file_path, 'r') as file:
        while True:
            # Read the number of atoms
            num_atoms_line = file.readline()
            if not num_atoms_line:
                break  # End of file
            num_atoms = int(num_atoms_line.strip())

            # Read the comment line
            comment_line = file.readline()
            if not comment_line:
                break  # End of file

            # Read the atomic data lines
            atoms = []
            for _ in range(num_atoms):
                line = file.readline()
                if not line:
                    break  # End of file or corrupted file
                atoms.append(line.strip())

            frames.append((num_atoms, comment_line.strip(), atoms))
    
    return frames

def write_xyz(file_path, frames):
    with open(file_path, 'w') as file:
        for num_atoms, comment, atoms in frames:
            file.write(f"{num_atoms}\n")
            file.write(f"{comment}\n")
            for atom in atoms:
                file.write(f"{atom}\n")

def filter_frames(frames, max_atoms):
    return [frame for frame in frames if frame[0] <= max_atoms]

def main(input_file, output_file, max_atoms=150):
    frames = read_xyz(input_file)
    filtered_frames = filter_frames(frames, max_atoms)
    write_xyz(output_file, filtered_frames)

if __name__ == "__main__":
    input_file = "input.xyz"
    output_file = "output.xyz"
    main(input_file, output_file)

