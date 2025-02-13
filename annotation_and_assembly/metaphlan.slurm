#!/bin/bash
#SBATCH -J metaphlan
#SBATCH -o /dm_data/zhaixp/20241107-dyn-soil/metaphlan.out
#SBATCH -e /dm_data/zhaixp/20241107-dyn-soil/metaphlan.err
#SBATCH -p dm_pub_cpu
#SBATCH -N 1
#SBATCH -n 32
#SBATCH --time=3-00:00:00


# 输入目录和输出目录
input_dir="/dm_data/zhaixp/20241107-dyn-soil/bbmap_output_many_seqs"
output_dir="/dm_data/zhaixp/20241107-dyn-soil/metaphlan_output/"

# 创建输出目录如果它不存在
mkdir -p "$output_dir"

# Bowtie2数据库和索引文件路径
bowtie2db="/dm_data/zhaixp/database/metaphlan_db"
index="mpa_vJun23_CHOCOPhlAnSGB_202403"

# 遍历输入目录下的所有.fasta文件
for fa_file in "$input_dir"/*.fastq; do
    # 提取文件名（不包括路径和扩展名）
    filename=$(basename "$fa_file" .fastq)
    
    # 定义输出文件路径
    output_file="$output_dir/${filename}_metaphlan.txt"
    
    # 执行metaphlan命令
    metaphlan --input_type fastq "$fa_file" \
        --bowtie2db "$bowtie2db" \
        --index "$index" \
        --offline \
        --nproc 32 \
        --output "$output_file"
done
