# 文件名: vasp2dp.py
# 运行方法: python vasp2dp.py
# 功能描述: 从多个系统中读取OUTCAR文件并将其转换为DeepMD格式数据

from dpdata import LabeledSystem, MultiSystems
from glob import glob

def convert_to_deepmd(pattern, output_dir, set_size=25):
    """
    处理多个系统并将其数据转换为DeepMD格式。

    参数:
    pattern (str): 匹配OUTCAR文件的全局模式。
    output_dir (str): 输出目录，用于存储DeepMD格式数据。
    set_size (int): DeepMD数据集的大小。
    """
    # 匹配所有符合模式的OUTCAR文件
    file_paths = glob(pattern)

    # 初始化MultiSystems对象
    multi_systems = MultiSystems()

    # 处理每个OUTCAR文件
    for file_path in file_paths:
        try:
            # 创建LabeledSystem对象
            labeled_system = LabeledSystem(file_path)
        except Exception as e:
            # 如果读取失败，输出错误信息
            print(f"Error processing file {file_path}: {e}")
            continue

        # 如果LabeledSystem对象不为空，添加到MultiSystems对象中
        if len(labeled_system) > 0:
            multi_systems.append(labeled_system)

    # 将MultiSystems对象转换为DeepMD的raw格式
    multi_systems.to_deepmd_raw(output_dir)

    # 将MultiSystems对象转换为DeepMD的npy格式
    multi_systems.to_deepmd_npy(output_dir, set_size=set_size)

    print(f"Processing complete. Data saved in '{output_dir}'.")

if __name__ == '__main__':
    # 全局模式，匹配OUTCAR文件的路径
    file_pattern = './*/*/*/OUTCAR'  # 请根据实际情况修改此路径

    # 输出目录
    output_directory = 'deepmd'

    # DeepMD数据集大小
    dataset_size = 25

    # 调用函数处理多系统数据
    convert_to_deepmd(file_pattern, output_directory, dataset_size)
