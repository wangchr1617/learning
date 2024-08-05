# 文件名: calculate_rmse.py
# 运行方法: python calculate_rmse.py
# 功能描述: 计算训练集中的能量、力和应力的均方根误差（RMSE）。

import numpy as np

def calculate_rmse_energy(file_path):
    """
    计算能量的均方根误差（RMSE）。

    参数:
        file_path (str): 能量数据文件路径。

    返回:
        float: 能量的RMSE，单位为meV/atom。
    """
    try:
        energy_data = np.loadtxt(file_path)
        rmse_energy = np.sqrt(np.mean((energy_data[:, 0] - energy_data[:, 1]) ** 2))
        return rmse_energy * 1000  # 单位转换为meV/atom
    except IOError:
        print(f"Error: Could not read file {file_path}")
    except ValueError:
        print("Error: File format is not correct")
    return None

def calculate_rmse_force(file_path):
    """
    计算力的均方根误差（RMSE）。

    参数:
        file_path (str): 力数据文件路径。

    返回:
        float: 力的RMSE，单位为meV/Å。
    """
    try:
        force_data = np.loadtxt(file_path)
        rmse_force = np.sqrt(np.mean((force_data[:, 3:6] - force_data[:, 0:3]) ** 2))
        return rmse_force * 1000  # 单位转换为meV/Å
    except IOError:
        print(f"Error: Could not read file {file_path}")
    except ValueError:
        print("Error: File format is not correct")
    return None

def calculate_rmse_virial(file_path):
    """
    计算应力的均方根误差（RMSE）。

    参数:
        file_path (str): 应力数据文件路径。

    返回:
        float: 应力的RMSE，单位为meV/atom。
    """
    try:
        virial_data = np.loadtxt(file_path)
        rmse_virial = np.sqrt(np.mean((virial_data[:, 6:12] - virial_data[:, 0:6]) ** 2))
        return rmse_virial * 1000  # 单位转换为meV/atom
    except IOError:
        print(f"Error: Could not read file {file_path}")
    except ValueError:
        print("Error: File format is not correct")
    return None

def main():
    """
    主函数，计算并打印能量、力和应力的RMSE。
    """
    # 计算能量的RMSE
    rmse_energy = calculate_rmse_energy('energy_train.out')
    if rmse_energy is not None:
        print(f"Energy RMSE: {rmse_energy:.2f} meV/atom")

    # 计算力的RMSE
    rmse_force = calculate_rmse_force('force_train.out')
    if rmse_force is not None:
        print(f"Force RMSE: {rmse_force:.2f} meV/Å")

    # 计算应力的RMSE
    rmse_virial = calculate_rmse_virial('virial_train.out')
    if rmse_virial is not None:
        print(f"Virial RMSE: {rmse_virial:.2f} meV/atom")

if __name__ == '__main__':
    main()
