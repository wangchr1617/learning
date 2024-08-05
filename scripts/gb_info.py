# 文件名: gb_info.py
# 运行方法: python gb_info.py 1 1 0 30
# 功能描述: 使用aimsgb库获取给定晶界方向和最大西格玛值的晶界信息。

import sys
from aimsgb import GBInformation

def get_gb_info(direction, sigma_max):
    """
    获取晶界信息。

    参数:
        direction (list): 晶界方向向量。
        sigma_max (int): 最大西格玛值。

    返回:
        GBInformation: 包含晶界信息的对象。
    """
    gb_info = GBInformation(direction, sigma_max, specific=False)
    return gb_info

def main():
    """
    主函数，解析命令行参数并输出晶界信息。
    """
    if len(sys.argv) != 5:
        print("Usage: python gb_info.py <x1> <x2> <x3> <sigma_max>")
        sys.exit(1)

    # 获取命令行参数
    x1 = int(sys.argv[1])
    x2 = int(sys.argv[2])
    x3 = int(sys.argv[3])
    sigma_max = int(sys.argv[4])

    # 获取晶界信息
    gb_info = get_gb_info([x1, x2, x3], sigma_max)

    # 输出晶界信息
    print(gb_info)
    print("Axis:", gb_info.axis)

if __name__ == '__main__':
    main()
