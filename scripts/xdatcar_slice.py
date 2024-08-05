"""
用法: python xdatcar_slice.py [begin] [end] [step]

功能描述:
该脚本从VASP输出的XDATCAR文件中读取特定范围的轨迹，并将其写入新的XDATCAR文件。

参数说明:
- [begin]: 轨迹的起始索引（包含）。
- [end]: 轨迹的结束索引（不包含）。
- [step]: 可选，轨迹的步长。

结果:
生成的轨迹将保存为新的XDATCAR文件，文件名为XDATCAR_REV。
"""

from ase.io import read, write
import sys

def main():
    # 获取命令行参数
    begin = int(sys.argv[1])
    end = int(sys.argv[2])
    
    # 尝试获取步长参数，如果没有提供则默认为None
    try:
        step = int(sys.argv[3])
    except IndexError:
        step = None

    # 从XDATCAR文件中读取指定范围的轨迹
    traj = read('./XDATCAR', index=slice(begin, end, step))
    print('The length of traj is:', len(traj))

    # 将轨迹写入新的XDATCAR文件
    write('./XDATCAR_REV', traj, format='vasp-xdatcar')
    print('New XDATCAR_REV file created successfully.')

if __name__ == "__main__":
    main()
