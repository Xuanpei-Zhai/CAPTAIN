#!/bin/bash
#SBATCH -J step1_basecalling
#SBATCH -o /dm_data/zhaixp/step1_basecalling_out1
#SBATCH -e /dm_data/zhaixp/step1_basecalling_error1
#SBATCH -p dm_pub_gpu
#SBATCH -N 1 
#SBATCH -n 17
#SBATCH --time=3-00:00:00
#SBATCH --gres=gpu:1 


/home_data/home/spst/zhaixp2022/dorado-0.8.3-linux-x64/bin/dorado basecaller /home_data/home/spst/zhaixp2022/dorado-0.8.3-linux-x64/bin/dna_r10.4.1_e8.2_400bps_sup@v5.0.0/ /dm_data/zhaixp/20250122-dyn/folder2/ --emit-fastq > /dm_data/zhaixp/20250122-dyn/calls3.fastq
