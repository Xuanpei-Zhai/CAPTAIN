#!/bin/bash
#SBATCH -J step2_filt
#SBATCH -o /dm_data/zhaixp/step2_filt_out
#SBATCH -e /dm_data/zhaixp/step2_filt_error
#SBATCH -p dm_pub_cpu
#SBATCH -N 1
#SBATCH -n 32
#SBATCH --time=3-00:00:00
#SBATCH --exclude=dm-c02n03

source ~/CAPTAIN/workflow/config.sh
source ~/anaconda3/etc/profile.d/conda.sh
conda activate captain

cd "$WORK_DIR" || exit 1

# 过滤和拆分FASTQ文件
NanoFilt -q 15 -l 20 --maxlength 3000 calls.fastq > calls_filt.fastq
seqkit split2 -f --quiet -j 16 -O fastq -s 10000 calls_filt.fastq

# 创建目标文件夹
for barcode in {01..10}; do
  for cpu in {1..10}; do
    mkdir -p "$WORK_DIR/barcode$barcode/fastq/cpu$cpu"
  done
done

# 获取所有FASTQ文件
SOURCE_DIR="$WORK_DIR/fastq"
files=($(ls $SOURCE_DIR/*.fastq))
total_files=${#files[@]}
files_per_barcode=$((total_files / 10))
files_per_cpu=$((files_per_barcode / 10))

# 分配文件
for (( i=0; i<total_files; i++ )); do
  barcode_index=$((i / files_per_barcode + 1))
  cpu_index=$(( (i % files_per_barcode) / files_per_cpu + 1 ))
  
  if (( barcode_index > 10 )); then
    barcode_index=10
  fi
  if (( cpu_index > 10 )); then
    cpu_index=10
  fi

  mv "${files[$i]}" "$WORK_DIR/barcode$(printf "%02d" $barcode_index)/fastq/cpu$cpu_index/"
done

echo "Files distribution completed."
