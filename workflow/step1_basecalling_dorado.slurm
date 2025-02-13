#!/bin/bash
#SBATCH -J step1_basecalling
#SBATCH -o /dm_data/zhaixp/step1_basecalling_out1
#SBATCH -e /dm_data/zhaixp/step1_basecalling_error1
#SBATCH -p dm_pub_gpu
#SBATCH -N 1 
#SBATCH -n 17
#SBATCH --time=3-00:00:00
#SBATCH --gres=gpu:1 

source ~/CAPTAIN/workflow/config.sh

$DORADO_PATH basecaller $KIT_PATH $POD5_DIR --emit-fastq > "$WORK_DIR/calls.fastq"
