def remove_config_frames(file_name):
    with open(file_name, 'r') as file:
        lines = file.readlines()
    
    filtered_lines = []
    i = 0
    while i < len(lines):
        num_atoms = int(lines[i].strip())
        comment = lines[i + 1].strip()
        frame_end = i + 2 + num_atoms
        
        if 'config' not in comment:
            filtered_lines.append(lines[i])
            filtered_lines.append(lines[i + 1])
            filtered_lines.extend(lines[i + 2:frame_end])
        
        i = frame_end
    
    with open('filtered_' + file_name, 'w') as file:
        file.writelines(filtered_lines)

remove_config_frames('selected.xyz')

