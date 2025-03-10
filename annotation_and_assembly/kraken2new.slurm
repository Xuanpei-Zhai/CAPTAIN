#!/bin/bash
#SBATCH -J kraken2_array
#SBATCH -o /dm_data/zhaixp/20241129-dyn-ocean/log/kraken2_out_%A_%a.out
#SBATCH -e /dm_data/zhaixp/20241129-dyn-ocean/log/kraken2_error_%A_%a.err
#SBATCH -p dm_pub_cpu
#SBATCH -N 1
#SBATCH -n 32
#SBATCH --time=3-00:00:00
#SBATCH --array=1-200

# 设置工作目录
cd /dm_data/zhaixp/20241129-dyn-ocean/

# 定义输入和输出目录
dir_bbmap="./bbmap_output_many_seqs/"
dir_bracken="./bracken_output_many_seqs/"
dir_kraken2="./kraken2_output_many_seqs/"

# 创建输出文件夹
mkdir -p "$dir_kraken2"
mkdir -p "$dir_bracken"

# 生成未处理的文件列表
unprocessed_list="unprocessed_files.txt"
> "$unprocessed_list"  # 清空文件

# 遍历所有bbmap中的FastQ文件，检查对应的Bracken输出是否存在
for fastq in "${dir_bbmap}"*.fastq; do
    base=$(basename "$fastq" .fastq)
    bracken_file="${dir_bracken}${base}_bracken.txt"
    if [[ ! -f "$bracken_file" ]]; then
        echo "$fastq" >> "$unprocessed_list"
    fi
done

# 如果没有需要处理的文件则退出
if [[ ! -s "$unprocessed_list" ]]; then
    echo "All files have been processed. Exiting."
    exit 0
fi

# 计算当前任务需要处理的文件范围
start_index=$(( (SLURM_ARRAY_TASK_ID - 1) * 100 + 1 ))
end_index=$(( SLURM_ARRAY_TASK_ID * 100 ))

# 获取文件列表中的对应文件
sed -n "${start_index},${end_index}p" "$unprocessed_list" | while read -r fastq_file; do
    base=$(basename "$fastq_file" .fastq)
    echo "Processing $base"

    # 运行Kraken2
    kraken2 -db /dm_data/zhaixp/database/kraken2_db --threads 32 \
        --output "${dir_kraken2}${base}_kraken2_output.txt" \
        --report "${dir_kraken2}${base}_kraken2_report.txt" \
        "$fastq_file"
    
    # 运行Bracken
    bracken -d /dm_data/zhaixp/database/kraken2_db \
        -i "${dir_kraken2}${base}_kraken2_report.txt" \
        -o "${dir_bracken}${base}_bracken.txt" \
        -r 300 -l G -t 50
done