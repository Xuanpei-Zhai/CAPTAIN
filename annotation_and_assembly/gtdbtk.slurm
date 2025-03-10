#!/bin/bash
#SBATCH -J gtdbtk
#SBATCH -o /dm_data/zhaixp/20241112-lmz-ocean/gtdbtk.out
#SBATCH -e /dm_data/zhaixp/20241112-lmz-ocean/gtdbtk.err
#SBATCH -p dm_debug
#SBATCH -N 1
#SBATCH -n 32
#SBATCH --time=3-00:00:00

source ~/anaconda3/etc/profile.d/conda.sh
conda activate gtdbtk

# 输入目录和输出目录
input_dir="/dm_data/zhaixp/20241112-lmz-ocean/miniasm_fasta_output"
output_dir="/dm_data/zhaixp/20241112-lmz-ocean/gtdbtk_output/"


gtdbtk classify_wf --genome_dir "$input_dir" --out_dir "$output_dir" --extension fasta --cpus 32 --mash_db /dm_data/zhaixp/database/gtdbtk_db/release220