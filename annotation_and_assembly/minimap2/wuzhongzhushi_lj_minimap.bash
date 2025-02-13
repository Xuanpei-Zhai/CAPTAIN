#!/bin/bash 

dir="/dm_data/zhaixp/20241129-dyn-Archaea-70/bbmap_output_many_seqs" 

conda activate
cd /dm_data/zhaixp/20241129-dyn-Archaea-70/
mkdir minimap2_output_many_seqs

for file in "$dir"/*.fastq
do 
    base=$(basename $file)
    echo "Processing $base" 

    minimap2 -ax map-ont --secondary=no /home_data/home/spst/zhaixp2022/minimap2/mut_ab.fasta "$file" > "./minimap2_output_many_seqs/${base%.fastq}_minimap2_output.txt"

done

