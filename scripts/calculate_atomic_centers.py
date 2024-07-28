# Usage: python calculate_atomic_centers.py

def calculate_atomic_centers(filepath='POSCAR'):
    """
    读取POSCAR文件并计算每种元素类型的原子中心位置。

    参数:
        filepath (str): POSCAR文件的路径，默认值为'POSCAR'。

    输出:
        每种元素类型的原子中心位置。
    """
    with open(filepath, 'r') as f:
        content = f.readlines()

    # 获取元素名称和每种元素的原子数量
    elements = [i for i in content[5].split()] 
    N_type = [int(i) for i in content[6].split()]

    res = []
    index = 8
    for i in range(len(N_type)):
        center = [0, 0, 0]
        tmp = content[index:index + N_type[i]]
        index += N_type[i]
        tmp1 = []
        for j in tmp:
            tmp2 = [float(each) for each in j.split()]
            tmp1.append(tmp2)
            center[0] += tmp2[0]
            center[1] += tmp2[1]
            center[2] += tmp2[2]
        res.append(tmp1)
        center = [center[x] / N_type[i] for x in range(len(center))]
        print(elements[i], 'atom center:', center)

if __name__ == "__main__":
    calculate_atomic_centers('POSCAR')
