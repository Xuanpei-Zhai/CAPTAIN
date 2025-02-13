#!/bin/bash

source ~/CAPTAIN/workflow/config.sh

if [ "$MODE" == "drop" ]; then
    sbatch "$CODE_DIR/demux/demux_drop/demux1.slurm"
    sbatch "$CODE_DIR/demux/demux_drop/demux2.slurm"
    sbatch "$CODE_DIR/demux/demux_drop/demux3.slurm"
    sbatch "$CODE_DIR/demux/demux_drop/demux4.slurm"
    sbatch "$CODE_DIR/demux/demux_drop/demux5.slurm"
    sbatch "$CODE_DIR/demux/demux_drop/demux6.slurm"
    sbatch "$CODE_DIR/demux/demux_drop/demux7.slurm"
    sbatch "$CODE_DIR/demux/demux_drop/demux8.slurm"
    sbatch "$CODE_DIR/demux/demux_drop/demux9.slurm"
    sbatch "$CODE_DIR/demux/demux_drop/demux10.slurm"
elif [ "$MODE" == "split" ]; then
    sbatch "$CODE_DIR/demux/demux_split/demux1.slurm"
    sbatch "$CODE_DIR/demux/demux_split/demux2.slurm"
    sbatch "$CODE_DIR/demux/demux_split/demux3.slurm"
    sbatch "$CODE_DIR/demux/demux_split/demux4.slurm"
    sbatch "$CODE_DIR/demux/demux_split/demux5.slurm"
    sbatch "$CODE_DIR/demux/demux_split/demux6.slurm"
    sbatch "$CODE_DIR/demux/demux_split/demux7.slurm"
    sbatch "$CODE_DIR/demux/demux_split/demux8.slurm"
    sbatch "$CODE_DIR/demux/demux_split/demux9.slurm"
    sbatch "$CODE_DIR/demux/demux_split/demux10.slurm"
elif [ "$MODE" == "16s" ]; then
    sbatch "$CODE_DIR/demux/demux_16s/demux1.slurm"
    sbatch "$CODE_DIR/demux/demux_16s/demux2.slurm"
    sbatch "$CODE_DIR/demux/demux_16s/demux3.slurm"
    sbatch "$CODE_DIR/demux/demux_16s/demux4.slurm"
    sbatch "$CODE_DIR/demux/demux_16s/demux5.slurm"
    sbatch "$CODE_DIR/demux/demux_16s/demux6.slurm"
    sbatch "$CODE_DIR/demux/demux_16s/demux7.slurm"
    sbatch "$CODE_DIR/demux/demux_16s/demux8.slurm"
    sbatch "$CODE_DIR/demux/demux_16s/demux9.slurm"
    sbatch "$CODE_DIR/demux/demux_16s/demux10.slurm"
else
    echo "Unknown MODE: $MODE. Please set MODE to 'drop' or 'split'."
    exit 1
fi
