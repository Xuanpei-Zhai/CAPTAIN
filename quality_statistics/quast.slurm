#!/bin/bash
#SBATCH -J quast_batch
#SBATCH -o /dm_data/zhaixp/20241129-dyn-ocean/quast_out.log
#SBATCH -e /dm_data/zhaixp/20241129-dyn-ocean/quast_error.log
#SBATCH -p dm_pub_cpu
#SBATCH -N 1
#SBATCH -n 32
#SBATCH --time=3-00:00:00
#SBATCH --exclude=dm-c02n03

# 设置工作目录
cd /dm_data/zhaixp/20241129-dyn-ocean/

# 创建输出目录（如果不存在）
mkdir -p quast_output

# 遍历 miniasm_fasta_output 目录中的所有 .fasta 文件
for fasta_file in ./miniasm_fasta_output/*.fasta; do
    # 获取文件名（不带路径和扩展名）
    base_name=$(basename "$fasta_file" .fasta)

    # 设置对应的输出目录
    output_dir=./quast_output/"$base_name"

    # 设置 nanopore 对应的 fastq 文件路径
    nanopore_file=./bbmap_output_many_seqs/"$base_name".fastq

    # 检查 nanopore 文件是否存在
    if [ -f "$nanopore_file" ]; then
        echo "Processing $fasta_file with Nanopore file $nanopore_file..."
        
        # 运行 quast 命令
        quast "$fasta_file" -o "$output_dir" -t 32 -f --nanopore "$nanopore_file"
    else
        echo "Warning: Nanopore file $nanopore_file not found. Skipping $fasta_file."
    fi
done


