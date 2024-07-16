# Usage: for i in */; do cd $i; python filter_xyz.py; cat new_dump.xyz >> ../iter02.xyz; cd ..; done
def parse_lattice_line(line):
    parts = line.split('Lattice="')[1].split('"')[0].split()
    lattice_params = [float(x) for x in parts]
    return lattice_params

def parse_force_lines(lines):
    forces = []
    for line in lines:
        force_values = [float(x) for x in line.split()[4:7]]
        forces.append(force_values)
    return forces

def should_delete_frame(lattice_params, forces, force_threshold, lattice_threshold):
    if any(abs(param) > lattice_threshold for param in lattice_params):
        return True
    for force in forces:
        if any(abs(f) > force_threshold for f in force):
            return True
    return False

def process_xyz_file(input_filename, output_filename, force_threshold, lattice_threshold):
    with open(input_filename, 'r') as infile, open(output_filename, 'w') as outfile:
        while True:
            line = infile.readline()
            if not line:
                break 
            atom_count = int(line.strip())
            info_line = infile.readline()
            lattice_params = parse_lattice_line(info_line)
            atom_lines = [infile.readline() for _ in range(atom_count)]
            forces = parse_force_lines(atom_lines)
            if not should_delete_frame(lattice_params, forces, force_threshold, lattice_threshold):
                outfile.write(line)
                outfile.write(info_line)
                outfile.writelines(atom_lines)

def main():
    input_filename = 'dump.xyz'
    output_filename = 'new_dump.xyz'
    force_threshold = 100
    lattice_threshold = 25
    process_xyz_file(input_filename, output_filename, force_threshold, lattice_threshold)

if __name__ == "__main__":
    main()
 
