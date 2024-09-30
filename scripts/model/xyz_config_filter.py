# 文件名: xyz_config_filter.py
# 运行方法: python xyz_config_filter.py
# 功能描述: 从XYZ文件中移除包含特定注释（如"config"）的帧，并将过滤后的结果保存到新的XYZ文件中。

def remove_config_frames(file_name):
    """
    从XYZ文件中移除包含特定注释的帧。

    参数:
        file_name (str): 要处理的XYZ文件名。
    """
    # 打开文件并读取所有行
    with open(file_name, 'r') as file:
        lines = file.readlines()
    
    filtered_lines = []  # 存储过滤后的行
    i = 0
    while i < len(lines):
        num_atoms = int(lines[i].strip())  # 读取原子数
        comment = lines[i + 1].strip()  # 读取注释行
        frame_end = i + 2 + num_atoms  # 计算帧的结束行索引
        
        # 检查注释中是否包含'config'
        if 'config' not in comment:
            # 如果不包含，则将该帧的所有行加入过滤后的列表
            filtered_lines.append(lines[i])
            filtered_lines.append(lines[i + 1])
            filtered_lines.extend(lines[i + 2:frame_end])
        
        # 移动到下一个帧的起始行
        i = frame_end
    
    # 将过滤后的内容写入新的文件
    output_file_name = 'filtered_' + file_name
    with open(output_file_name, 'w') as file:
        file.writelines(filtered_lines)

    print(f"Filtered data saved to {output_file_name}")

# 调用函数处理文件
remove_config_frames('selected.xyz')
