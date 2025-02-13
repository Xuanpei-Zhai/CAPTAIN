#!/bin/bash
#SBATCH -J kraken2
#SBATCH -o /dm_data/zhaixp/20241217-dyn-soil/kraken2.out
#SBATCH -e /dm_data/zhaixp/20241217-dyn-soil/kraken2.err
#SBATCH -p dm_pub_cpu
#SBATCH -N 1
#SBATCH -n 32
#SBATCH --time=3-00:00:00


dir="/dm_data/zhaixp/20241217-dyn-soil/bbmap_output_many_seqs/" 


cd /dm_data/zhaixp/20241217-dyn-soil/ 
mkdir kraken2_output_many_seqs
mkdir bracken_output_many_seqs

for file in ${dir}*.fastq
do 
    base=$(basename $file)
    echo "Processing $base" 


    kraken2 -db /dm_data/zhaixp/database/kraken2_db --threads 32 --output "./kraken2_output_many_seqs/${base%.fastq}_kraken2_output.txt" --report "./kraken2_output_many_seqs/${base%.fastq}_kraken2_report.txt" "$file" 
    bracken -d /dm_data/zhaixp/database/kraken2_db -i "./kraken2_output_many_seqs/${base%.fastq}_kraken2_report.txt" -o "./bracken_output_many_seqs/${base%.fastq}_bracken.txt" -r 300 -l G -t 50
done


