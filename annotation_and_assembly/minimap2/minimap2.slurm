#!/bin/bash
#SBATCH -J minimap
#SBATCH -o /dm_data/zhaixp/20241129-dyn-Archaea-70/minimap_o
#SBATCH -e /dm_data/zhaixp/20241129-dyn-Archaea-70/minimap_e
#SBATCH -p dm_pub_cpu
#SBATCH -N 1
#SBATCH -n 32
#SBATCH --time=3-00:00:00
#SBATCH --exclude=dm-c02n03

/home_data/home/spst/zhaixp2022/CAPTAIN/annotation_and_assembly/minimap2/wuzhongzhushi_lj_minimap.bash
