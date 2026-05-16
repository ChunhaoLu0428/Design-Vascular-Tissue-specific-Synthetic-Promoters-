#!/bin/bash

# 定义基序文件路径 τ=100、099、098、097、096、095、090、085
motif_file="/path/to/each/τ/NRmerged_streme.meme"

# 定义输入的序列目录
input_base_dir="path/to/each/length/promoter_seq/"

# 定义基序文件的输出根目录
fimo_output_base_dir="/path/to/each/τ/fimo_out"

# 遍历所有子目录（遍历 LM50/ 下的所有子目录）
for sub_dir in "$input_base_dir"/*/; do
    # 确定当前子目录的名称，例如 "500" 或其他子目录名称
    sub_dir_name=$(basename "$sub_dir")
    
    # 检查每个子目录下的 promoter_output 是否存在 merged_promoters.fa
    promoter_file="$sub_dir/promoter_output/merged_promoters.fa"
    
    if [[ -f "$promoter_file" ]]; then
        # 设置 fimo 的输出目录，根据子目录名称决定
        output_dir="$fimo_output_base_dir/$sub_dir_name"
        
        # 创建输出目录（如果不存在）
        mkdir -p "$output_dir"
        
        # 执行 fimo 命令
        echo "Running fimo for $promoter_file in $sub_dir_name..."
        nohup fimo --oc "$output_dir" "$motif_file" "$promoter_file" &
    else
        echo "File $promoter_file does not exist. Skipping..."
    fi
done

# 等待所有后台任务完成
wait
echo "FIMO analysis completed for all directories."
