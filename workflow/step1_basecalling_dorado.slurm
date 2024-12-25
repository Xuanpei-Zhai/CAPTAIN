#!/bin/bash
#SBATCH -J step1_basecalling
#SBATCH -o /dm_data/zhaixp/20241217-dyn-soil/step1_basecalling_out
#SBATCH -e /dm_data/zhaixp/20241217-dyn-soil/step1_basecalling_error
#SBATCH -p dm_pub_gpu
#SBATCH -N 1 
#SBATCH -n 16
#SBATCH --time=3-00:00:00
#SBATCH --gres=gpu:1 

nvidia-smi

/home_data/home/spst/zhaixp2022/dorado-0.8.3-linux-x64/bin/dorado basecaller /home_data/home/spst/zhaixp2022/dorado-0.8.3-linux-x64/bin/dna_r10.4.1_e8.2_400bps_sup@v5.0.0/ /dm_data/zhaixp/20241217-dyn/soil-droplet-1th/20241217_1649_P2S-02431-B_PBA77234_52a4f57d/pod5 --emit-fastq > /dm_data/zhaixp/20241217-dyn-soil/calls.fastq
