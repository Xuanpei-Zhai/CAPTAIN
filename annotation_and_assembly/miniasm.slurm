#!/bin/bash
#SBATCH -J miniasm_array
#SBATCH -o /dm_data/zhaixp/20241217-dyn-soil/log/miniasm_out_%A_%a.out
#SBATCH -e /dm_data/zhaixp/20241217-dyn-soil/log/miniasm_error_%A_%a.err
#SBATCH -p dm_pub_cpu
#SBATCH -N 1
#SBATCH -n 32
#SBATCH --time=3-00:00:00
#SBATCH --exclude=dm-c02n03
#SBATCH --array=1-999

# 设置工作目录
cd /dm_data/zhaixp/20241217-dyn-soil/

# 创建输出文件夹
mkdir -p miniasm_paf_output
mkdir -p miniasm_fasta_output

# 计算当前任务需要处理的文件范围
start_index=$(( (SLURM_ARRAY_TASK_ID - 1) * 10 + 1 ))
end_index=$(( SLURM_ARRAY_TASK_ID * 10 ))

# 获取文件列表中的对应文件
sed -n "${start_index},${end_index}p" file_list.txt | while read -r fastq_file; do
    # 获取文件名（去掉路径和后缀）
    base_name=$(basename "$fastq_file" .fastq)

    # 执行minimap2命令，生成paf文件
    minimap2 -x ava-ont -t 32 "$fastq_file" "$fastq_file" | gzip -1 > ./miniasm_paf_output/"$base_name".paf.gz

    # 执行miniasm命令，生成gfa文件
    miniasm -f "$fastq_file" ./miniasm_paf_output/"$base_name".paf.gz -c 1 -m 100 -s 100 -g 500 -h 10000 -I 0.6 -R > ./miniasm_paf_output/"$base_name".gfa

    # 提取fasta序列并保存到fasta_output文件夹中
    awk '/^S/{print ">"$2"\n"$3}' ./miniasm_paf_output/"$base_name".gfa | fold > ./miniasm_fasta_output/"$base_name".fasta
    rm ./miniasm_paf_output/"$base_name".paf.gz
    rm ./miniasm_paf_output/"$base_name".gfa
done

