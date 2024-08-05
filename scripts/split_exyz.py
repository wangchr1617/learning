"""
用法: python split_exyz.py -i [input_file] --ranges [frame_ranges] [--nframes N] [--perrange]

功能描述:
该脚本用于分割EXYZ文件，选择指定范围或数量的帧，并输出选定的帧。

参数说明:
- -i, --input: 输入的EXYZ文件路径（必需）。
- --ranges: 帧的范围，格式为"0:2,4,7:9"，表示选择帧[0,1,4,7,8]（必需）。
- --nframes: 要随机选择的帧数。
- --perrange: 如果设置，将从每个范围中选择NFRAMES帧。

注意事项:
- 确保帧索引和范围格式正确。
- 参数NFRAMES应为正整数。
"""

import io as _io
import sys as _sys
import random as _rn
import argparse as _ar

def parse_args() -> _ar.Namespace:
    """
    解析命令行参数并对其进行验证。

    返回:
    argparse.Namespace: 包含解析后参数的命名空间。
    """
    description = 'split_exyz.py: Splitting an EXYZ file.'
    parser = _ar.ArgumentParser(description=description)
    parser.add_argument('-i', '--input', required=True, type=str,
                        help='The input EXYZ file.')
    parser.add_argument('--ranges', required=True, type=str,
                        help='The frame ranges. For example, "0:2,4,7:9" ' +
                             'will select frames: [0,1,4,7,8]')
    parser.add_argument('--nframes', type=int,
                        help='Number of frames to randomly select.')
    parser.add_argument('--perrange', action='store_true',
                        help='Select NFRAMES from each range.')
    
    if len(_sys.argv) == 1:
        parser.print_help()
        exit()
        
    return parser.parse_args()

def read_frame(fp: _io.TextIOWrapper, n_frame: int) -> tuple:
    """
    从输入流中读取一个帧。

    参数:
    fp (io.TextIOWrapper): 文件指针。
    n_frame (int): 当前帧的编号。

    返回:
    tuple: 包含原子数量和帧数据的元组。
    """
    lines = []
    n_atoms = fp.readline()
    
    if n_atoms == '':
        return (0, [])
    else:
        n_atoms = int(n_atoms)
        
    for i in range(n_atoms + 1):
        line = fp.readline()
        if i == 0:
            if not ('energy' in line or 'energies' in line):
                message = f'Frame {n_frame} contains a bad comment line: "{line.strip()}".'
                raise RuntimeError(message)
        else:
            words = line.strip().split(' ')
            if (len(words) - words.count('')) < 4:
                message = f'Frame {n_frame} contains a bad data line: "{line.strip()}". '
                message += 'This may be caused by a mismatch between the given atom number and the data line number.'
                raise RuntimeError(message)
        if line != '':
            lines.append(line.strip())
            
    if len(lines) != (n_atoms + 1):
        message = f'Mismatch between atom number {n_atoms} and xyz data:\n'
        for line in lines:
            message += line + '\n'
        raise RuntimeError(message)
        
    return (n_atoms, lines)

def parse_input_frame_list(frame_list: str) -> list:
    """
    解析帧选择字符串。

    参数:
    frame_list (str): 帧选择字符串。

    返回:
    list: 解析后的帧列表。
    """
    result = []
    
    for s in frame_list.split(','):
        if ':' not in s:
            if not s.isnumeric():
                raise SyntaxError('Unknown frame index: ' + s)
            result.append(int(s))
        else:
            l = s.split(':')
            if not l[0].isnumeric() or len(l) != 2 or not l[1].isnumeric():
                raise SyntaxError('Unknown frame range: ' + s)
            if int(l[0]) >= int(l[1]):
                message = f'Bad frame range [{int(l[0])}, {int(l[1])}]! The first index should be smaller than the second index!'
                raise RuntimeError(message)
            result.extend(list(range(int(l[0]), int(l[1]))))
    
    return result

def make_frame_list(frame_ranges: list, n_frames: int, per_range: bool) -> list:
    """
    获取要写入的帧的索引。

    参数:
    frame_ranges (list): 帧的范围列表。
    n_frames (int): 要随机选择的帧数。
    per_range (bool): 是否从每个范围中选择NFRAMES帧。

    返回:
    list: 要写入的帧的索引列表。
    """
    frame_list = []
    
    for element in frame_ranges:
        if isinstance(element, int):
            frame_list.append(element)
        elif isinstance(element, list):
            if per_range and n_frames is not None:
                if n_frames < len(element):
                    element = _rn.sample(element, n_frames)
            frame_list.extend(element)
            
    if n_frames is not None and not per_range:
        frame_list = _rn.sample(frame_list, n_frames)
        
    frame_list = list(set(frame_list))
    _rn.shuffle(frame_list)
    
    return frame_list

def print_frame(n_atoms: int, lines: list):
    """
    将一个帧打印到屏幕上。

    参数:
    n_atoms (int): 原子数量。
    lines (list): 帧数据的行列表。
    """
    print(n_atoms)
    for line in lines:
        print(line)

def loop(file_name: str, frame_list: list):
    """
    遍历输入文件。

    参数:
    file_name (str): 输入文件的名称。
    frame_list (list): 要输出的帧列表。
    """
    count = 0
    frames = []
    with open(file_name, 'r') as fp:
        while True:
            n_atoms, lines = read_frame(fp, count)
            if len(lines) == 0:
                break
            frames.append((n_atoms, lines))
            count += 1

    n_frames_total = len(frames)
    
    for i in frame_list:
        if i >= n_frames_total:
            message = f'Cannot read frame: {i}'
            raise RuntimeError(message)
        else:
            print_frame(frames[i][0], frames[i][1])
            
    del frames

if __name__ == '__main__':
    # 解析命令行参数
    parser = parse_args()
    
    # 验证参数nframes
    if parser.nframes is not None and parser.nframes <= 0:
        raise RuntimeError('The NFRAMES parameter should be positive!')
    
    # 解析和生成帧列表
    frame_list = parse_input_frame_list(parser.ranges)
    frame_list = make_frame_list(frame_list, parser.nframes, parser.perrange)
    
    # 输出选定的帧
    print("Selected frames:", file=_sys.stderr)
    for i in frame_list:
        print(i, file=_sys.stderr)
    
    # 处理输入文件并输出选定帧
    loop(parser.input, frame_list)
