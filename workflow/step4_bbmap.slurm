#!/bin/bash 
#SBATCH -J step4_demux
#SBATCH -o /dm_data/zhaixp/step4_demux_out
#SBATCH -e /dm_data/zhaixp/step4_demux_error
#SBATCH -p dm_debug
#SBATCH -N 1
#SBATCH -n 32
#SBATCH --time=1-00:00:00
#SBATCH --exclude=dm-c02n03 

source ~/CAPTAIN/workflow/config.sh
source ~/anaconda3/etc/profile.d/conda.sh
conda activate captain

cd "$WORK_DIR"

find ./barcode*/cpu*/barcode_match*_output/process_* -type f -name "*.fastq" -print0 | xargs -0 cat >> ./bbmap_input_nocut.fastq


cutadapt -g AGATGTGTATAAGAGACAG -O 9 -a CTGTCTCTTATACACATCT -O 10 -m 50 --cores 32 -e 0.15 -o ./bbmap_input.fastq ./bbmap_input_nocut.fastq 


# 设置输入和输出目录
INPUT_FILE="$WORK_DIR/bbmap_input.fastq"
OUTPUT_DIR="$WORK_DIR/bbmap_output"
SPLIT_DIR="$WORK_DIR/split_fastq"

# 创建拆分文件的目录
mkdir -p $SPLIT_DIR
mkdir -p $OUTPUT_DIR

# 确定拆分文件的数量
NUM_FILES=30

# 获取文件中的总行数
TOTAL_LINES=$(wc -l < $INPUT_FILE)

# 因为FASTQ格式每个记录是4行，所以每个拆分文件应该有大致相同的记录数
LINES_PER_FILE=$((TOTAL_LINES / NUM_FILES / 4 * 4))

# 拆分文件
split -l $LINES_PER_FILE -d $INPUT_FILE $SPLIT_DIR/split_

# 对每个拆分后的文件运行demuxbyname
for FILE in $SPLIT_DIR/split_*
do
    demuxbyname.sh in=$FILE out=$OUTPUT_DIR/out_%.fastq prefixmode=f delimiter=_ -Xmx30g ow=f
done


python $CODE_DIR/demux/statistical_bbmap.py -i $WORK_DIR/bbmap_output -o $WORK_DIR/ -s 1000
