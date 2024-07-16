import sys
from ase.io import read
from ase.geometry import get_angles

def calculate_angle(poscar_path, atom_indices):
    # 加载POSCAR文件
    atoms = read(poscar_path)

    # 确保输入的原子索引是有效的
    if all(i < len(atoms) for i in atom_indices):
        # 计算并返回键角
        a = atom_indices[0]
        b = atom_indices[1]
        c = atom_indices[2] 
        v1 = (atoms.positions[a] - atoms.positions[b]).reshape(-1,3)
        v2 = (atoms.positions[c] - atoms.positions[b]).reshape(-1,3)
        angle = get_angles(v1, v2, cell=None, pbc=None)
        return angle
    else:
        raise IndexError("指定的原子索引超出了范围。")

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("使用方法: python script.py <POSCAR路径> <原子索引1> <原子索引2> <原子索引3>")
    else:
        poscar_path = sys.argv[1]
        atom_indices = [int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])]
        angle = calculate_angle(poscar_path, atom_indices)
        print(angle[0])
