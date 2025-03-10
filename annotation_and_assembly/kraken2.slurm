#!/bin/bash
#SBATCH -J kraken2
#SBATCH -o /dm_data/zhaixp/20241203-dyn-ocean/kraken2.out
#SBATCH -e /dm_data/zhaixp/20241203-dyn-ocean/kraken2.err
#SBATCH -p dm_pub_cpu
#SBATCH -N 1
#SBATCH -n 32
#SBATCH --time=3-00:00:00

dir_bbmap="/dm_data/zhaixp/20241203-dyn-ocean/bbmap_output_many_seqs/"
dir_bracken="/dm_data/zhaixp/20241203-dyn-ocean/bracken_output_many_seqs/"

cd /dm_data/zhaixp/20241203-dyn-ocean/ 
mkdir -p kraken2_output_many_seqs
mkdir -p bracken_output_many_seqs

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

echo "Processing files listed in $unprocessed_list"

# 根据列表处理未处理的FastQ文件
while IFS= read -r file; do
    base=$(basename "$file" .fastq)
    echo "Processing $base"

    # 运行Kraken2
    kraken2 -db /dm_data/zhaixp/database/kraken2_db --threads 32 \
        --output "./kraken2_output_many_seqs/${base}_kraken2_output.txt" \
        --report "./kraken2_output_many_seqs/${base}_kraken2_report.txt" \
        --memory-mapping "$file"
    
    # 运行Bracken
    bracken -d /dm_data/zhaixp/database/kraken2_db \
        -i "./kraken2_output_many_seqs/${base}_kraken2_report.txt" \
        -o "./bracken_output_many_seqs/${base}_bracken.txt" \
        -r 300 -l G -t 50
done < "$unprocessed_list"