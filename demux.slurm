#!/bin/bash 
#SBATCH -J step4_demux
#SBATCH -o /dm_data/zhaixp/20241129-dyn-Archaea-70/step4_demux_out
#SBATCH -e /dm_data/zhaixp/20241129-dyn-Archaea-70/step4_demux_error
#SBATCH -p dm_pub_cpu
#SBATCH -N 1
#SBATCH -n 32
#SBATCH --time=1-00:00:00
#SBATCH --exclude=dm-c02n03 

cd /dm_data/zhaixp/20241129-dyn-Archaea-70


# 设置输入和输出目录
INPUT_FILE="/dm_data/zhaixp/20241129-dyn-Archaea-70/bbmap_input.fastq"
OUTPUT_DIR="/dm_data/zhaixp/20241129-dyn-Archaea-70/bbmap_output"
SPLIT_DIR="/dm_data/zhaixp/20241129-dyn-Archaea-70/split_fastq"

# 对每个拆分后的文件运行demuxbyname
for FILE in $SPLIT_DIR/split_*
do
    demuxbyname.sh in=$FILE out=$OUTPUT_DIR/out_%.fastq prefixmode=f delimiter=_ -Xmx30g ow=f
done


python /home_data/home/spst/zhaixp2022/code/statistical_bbmap.py -i /dm_data/zhaixp/20241129-dyn-Archaea-70/bbmap_output -o /dm_data/zhaixp/20241129-dyn-Archaea-70/ -s 10000
