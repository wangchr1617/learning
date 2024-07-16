import os

def delete_files_containing_sb(directory):
    for filename in os.listdir(directory):
        if filename.endswith(".vasp"):
            filepath = os.path.join(directory, filename)
            try:
                with open(filepath, 'r', encoding='utf-8') as file:
                    contents = file.read()
                if 'Sb' in contents:
                    os.remove(filepath) 
                    print(f"Deleted '{filename}' as it contains 'Sb'.")
            except Exception as e:
                print(f"Failed to read or delete '{filename}'. Error: {e}")

current_directory = os.getcwd()
delete_files_containing_sb(current_directory)

