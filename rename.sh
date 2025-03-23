#!/bin/bash

# 遍历 1 到 25
for i in $(seq 1 25); do
    # 格式化文件名
    input_file="phase0.A.counts_AAGCTT.25g${i}.txt"
    output_file="${input_file}.1"

    # 运行 awk 命令
    awk '
        NR==FNR {a[$1]; next}  # 处理 common.name 文件，将第一列的字段存储在数组 a 中
        {
            # 检查第一列是否在数组 a 中
            if ($1 in a) {
                print
            }
        }
    ' 1.name "$input_file" > "$output_file"

    echo "Processed $input_file"
done

