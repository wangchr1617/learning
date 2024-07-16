# Usage: ./xyz_grouped.sh
sed -i 's/Properties=species:S:1:pos:R:3/Properties=species:S:1:pos:R:3:group:I:1/' model.xyz

input_file="model.xyz"
output_file="grouped_model.xyz"

# 读取输入文件，并处理行
while IFS= read -r line
do
    if [[ $line == Ge* ]]; then
        echo "$line  0" >> $output_file
    elif [[ $line == Te* ]]; then
        echo "$line  1" >> $output_file
    else
        echo "$line" >> $output_file
    fi
done < "$input_file"

mv $output_file $input_file

