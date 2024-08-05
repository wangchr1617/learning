# 文件名: gpumd_tool.py
# 运行方法: python gpumd_tool.py [job_type] [path]
# 功能描述: 使用GPUMD进行分子动力学模拟、主动学习和预测。

import matplotlib
import numpy as np
from scipy.spatial.distance import cdist
from pathlib import Path
from calorine.gpumd import *
from calorine.nep import get_descriptors
from ase.io import read as ase_read
from ase.io import write as ase_write
import matplotlib.pyplot as plt
from monty.os import cd
from sklearn.decomposition import PCA
import argparse
import datetime
import glob
import logging
import os
import sys
import shutil
import subprocess

# 设置matplotlib后端为Agg，以支持无界面图形保存
matplotlib.use("Agg")

# 配置日志输出
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S',
    stream=sys.stdout
)

# 采样间隔，每隔NumSamples抽取一个数据点
NumSamples = 2

def select(new_data, now_data=[], min_distance=None, min_select=1, max_select=None):
    """
    选择那些与给定数据最远的数据点。
    
    参数:
        new_data (2d list or array): 需要选择的数据点集。
        now_data (2d list or array): 已经存在的数据点集。默认为空。
        min_distance (float, optional): 
            如果两个点之间的距离超过最小距离，停止选择。默认为None。
        min_select (int, optional): 要选择的最小数据点数，这可能导致某些点之间的距离小于给定的最小距离。默认为1。
        max_select (int, optional): 要选择的最大数据点数。默认为None（无限制）。

    返回:
        list of int: 选择的数据点的索引列表。
    """
    metric = 'euclidean'
    metric_para = {}
    max_select = max_select or len(new_data)
    to_add = []

    if len(new_data) == 0:
        return to_add

    if len(now_data) == 0:
        to_add.append(0)
        now_data.append(new_data[0])

    # 计算新数据点与现有数据点之间的距离
    distances = np.min(cdist(new_data, now_data, metric=metric, **metric_para), axis=1)

    while np.max(distances) > min_distance or len(to_add) < min_select:
        i = np.argmax(distances)
        to_add.append(i)
        if len(to_add) >= max_select:
            break
        # 更新与新选择数据点的最小距离
        distances = np.minimum(distances, cdist([new_data[i]], new_data, metric=metric)[0])

    return to_add

def run(run_cmd: str, run_dir: Path):
    """
    执行指定的命令，并将输出和错误重定向到文件。

    参数:
        run_cmd (str): 要执行的命令。
        run_dir (Path): 执行命令的工作目录。
    """
    start = datetime.datetime.now()
    logging.info("\t开始计算")

    vasp_cmd = [os.path.expanduser(os.path.expandvars(run_cmd))]
    with cd(run_dir), open(f"{run_cmd}.out", "w") as f_std, open(f"{run_cmd}.err", "w", buffering=1) as f_err:
        subprocess.check_call(vasp_cmd, stdout=f_std, stderr=f_err)
    logging.info("\t计算完成" + f"\t耗时：{datetime.datetime.now() - start}")

def remove_garbage_structure(atoms_list):
    """
    删除崩溃的结构。
    
    参数:
        atoms_list (list): 包含Atoms对象的列表。

    返回:
        list: 过滤后的Atoms对象列表。
    """
    result = []
    for atoms in atoms_list:
        position = atoms.get_all_distances()
        # 检查最小距离，如果小于1，则认为是崩溃结构
        if (np.min(position[position > 0])) < 1:
            continue
        result.append(atoms)

    return result

def verify_path(path: Path) -> None:
    """
    验证路径是否存在，如果不存在则创建路径（支持多级目录）。

    参数:
        path (Path): 要验证或创建的路径。
    """
    if not path.exists():
        os.makedirs(path)

def cp_file(source_file: Path, destination_dir: Path) -> None:
    """
    复制文件到目标目录。

    参数:
        source_file (Path): 要复制的文件路径。
        destination_dir (Path): 目标目录路径。
    """
    src_files = glob.glob(source_file.as_posix())
    for i in src_files:
        logging.debug(f"\t复制文件：{i} -> {destination_dir.as_posix()}")
        shutil.copy(i, destination_dir.as_posix())

def iter_path(glob_strs: list):
    """
    装饰器，用于遍历路径并对文件进行处理。

    参数:
        glob_strs (list): 用于匹配文件的glob字符串列表。
    """
    def decorator(func):
        def wrapper(path: Path | str, *args, **kwargs):
            if isinstance(path, str):
                path = Path(path)
            if path.is_dir():
                parent = path
            else:
                parent = path.parent
            result =[]
            for glob_str in glob_strs:

                for i in parent.glob(glob_str):
                    if path.is_file():
                        if i.name != path.name:
                            continue

                    try:
                        result.append(func(i, *args, **kwargs))
                    except KeyboardInterrupt:
                        return
                    except Exception as e:
                        logging.error(e)
                        pass
            return result
        return wrapper

    return decorator

@iter_path(["*.xyz", "*.vasp"])
def molecular_dynamics(path: Path, temperature, run_time):
    """
    根据指定的文件夹，计算文件夹下所有的xyz或vasp文件。

    参数:
        path (Path): 输入文件路径。
        temperature (float): 模拟的温度。
        run_time (int): 运行时间。

    返回:
        Path: 计算结果的路径。
    """
    # 读取原子结构
    if path.suffix == ".vasp":
        atoms = ase_read(path, 0, format="vasp")
    else:
        atoms = ase_read(path, 0, format="extxyz")

    # 创建计算路径
    md_path = root_path.joinpath(f"cache/{atoms.symbols}/{run_time}/md-{temperature}k")
    verify_path(md_path)
    logging.info(f"路径：{md_path.as_posix()}")

    # 编写运行文件
    run_in = [('potential', 'nep.txt'),
              ('velocity', temperature),
              ('ensemble', ('nvt_nhc', temperature, temperature, '100')),
              ('time_step', 1.0),
              ('dump_thermo', 1000),
              ('dump_exyz', ('1000', '0', '0')),
              ('run', 1000 * run_time)]

    write_runfile(md_path.joinpath("run.in"), run_in)
    cp_file(root_path.joinpath("nep.txt"), md_path.joinpath("nep.txt"))
    atoms.write(md_path.joinpath("model.xyz"), format="extxyz")

    # 运行GPUMD
    run("gpumd", md_path)

    # 读取结果
    data = read_thermo(md_path.joinpath("thermo.out").as_posix(), len(atoms))
    potential_energy = data.potential_energy.to_numpy(dtype='float')

    # 绘制能量图
    fig = plt.figure()
    plt.plot(list(range(potential_energy.shape[0])), potential_energy)
    plt.rcParams['font.size'] = 14
    plt.rcParams['axes.unicode_minus'] = False
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'
    plt.savefig(md_path.joinpath("md_energy.png"), dpi=150)

    return md_path

def select_structures(train, new: Path, max_selected=20):
    """
    选择结构，去除崩溃结构并进行特征选择。

    参数:
        train (list): 训练数据的Atoms对象列表。
        new (Path): 新的Atoms对象文件路径。
        max_selected (int): 最大选择数目。

    返回:
        list: 选择的Atoms对象。
    """
    # 移除崩溃结构
    new_atoms = ase_read(new, ":", format="extxyz")
    new_atoms = remove_garbage_structure(new_atoms)

    # 获取特征描述符
    train_des = np.array([np.mean(get_descriptors(i, "nep.txt"), axis=0) for i in train])
    new_des = np.array([np.mean(get_descriptors(i, "nep.txt"), axis=0) for i in new_atoms])

    # 选择新结构
    selected_i = select(np.vstack([train_des, new_des]), train_des, min_distance=0.01, max_select=max_selected, min_select=0)

    # 绘制PCA图
    reducer = PCA(n_components=2)
    reducer.fit(new_des)
    proj = reducer.transform(new_des)
    fig = plt.figure()
    plt.scatter(proj[:, 0], proj[:, 1], label='all data')

    if selected_i:
        selected_proj = reducer.transform(np.array([new_des[i - train_des.shape[0]] for i in selected_i]))
        plt.scatter(selected_proj[:, 0], selected_proj[:, 1], label='selected data')
    plt.legend()
    plt.axis('off')
    plt.savefig(new.with_name('select.png'))

    return [new_atoms[i - train_des.shape[0]] for i in selected_i]

def auto_learn():
    """
    主动学习迭代。
    需要nep.txt、nep.in和train.xyz文件。
    """
    # 定义迭代时间，单位ps
    times = [10000]
    temperatures = range(50, 1000, 50)
    trainxyz = ase_read("train.xyz", ":", format="extxyz")
    
    for epoch, run_time in enumerate(times):
        logging.info(f"开始第{epoch + 1}次主动学习，采样时长：{run_time} ps。")
        
        # 存放每次epoch新增的训练集
        new_atoms = []

        # 进行gpumd采样
        for temperature in temperatures:
            logging.info(f"GPUMD采样中，温度：{temperature}k。时长：{run_time}ps")

            md_paths = molecular_dynamics("./s/", temperature=temperature, run_time=run_time)
            
            # 筛选出结构
            for md_path in md_paths:
                selected = select_structures(trainxyz, md_path.joinpath("dump.xyz"), max_selected=20)
                logging.info(f"得到{len(selected)}个结构")
                for i, atom in enumerate(selected):
                    atom.info["Config_type"] = f"epoch-{epoch + 1}-{run_time}ps-{temperature}k-{i + 1}"
                new_atoms.extend(selected)

        logging.info(f"本次主动学习新增了{len(new_atoms)}个结构。")

        # 保存新结构
        ase_write(root_path.joinpath(f"result/learn-epoch-{epoch}-{run_time}ps.xyz"), new_atoms, format="extxyz")
        break

        # 然后进行nep训练

def prediction(path: Path):
    """
    执行预测任务。

    参数:
        path (Path): 输入文件路径。
    """
    # 实现预测功能
    pass

def build_argparse():
    """
    构建命令行参数解析器。

    返回:
        ArgumentParser: 参数解析器。
    """
    parser = argparse.ArgumentParser(
        description="""GPUMD工具. 
        可以批量md和主动学习 """,
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_subparsers()

    parser.add_argument(
        "job_type", choices=["prediction", "md", "learn"], help="任务类型：预测、分子动力学模拟、主动学习"
    )
    parser.add_argument(
        "path", type=Path, help="要计算的xyz路径，或者要批量计算的文件夹。"
    )

    return parser

if __name__ == '__main__':
    # 采样
    parser = build_argparse()
    args = parser.parse_args()

    if not os.path.exists("./result"):
        os.mkdir("./result")
    root_path = Path("./")

    if args.job_type == "md":
        for t in range(50, 1000, 50):
            molecular_dynamics(args.path, temperature=t)
    elif args.job_type == "prediction":
        prediction(args.path)
    elif args.job_type == "learn":
        auto_learn()
