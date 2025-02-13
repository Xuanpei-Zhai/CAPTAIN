#!/bin/bash
#SBATCH -J step3_8_demux
#SBATCH -p dm_pub_cpu
#SBATCH -o /dm_data/zhaixp/step3_8_demux_out
#SBATCH -e /dm_data/zhaixp/step3_8_demux_error
#SBATCH -N 1
#SBATCH -n 32
#SBATCH --time=3-00:00:00
#SBATCH --exclude=dm-c02n03
source ~/CAPTAIN/workflow/config.sh
source ~/anaconda3/etc/profile.d/conda.sh
conda activate captain

python "$CODE_DIR"/search_barcode/search_fastq_dyn_pro.py -i "$WORK_DIR"/barcode01/fastq/cpu8 -o "$WORK_DIR"/barcode01/cpu8/search_output -s GTGATCTCGCTC -t 20
python "$CODE_DIR"/search_barcode/search_fastq_dyn_pro.py -i "$WORK_DIR"/barcode02/fastq/cpu8 -o "$WORK_DIR"/barcode02/cpu8/search_output -s GTGATCTCGCTC -t 20
python "$CODE_DIR"/search_barcode/search_fastq_dyn_pro.py -i "$WORK_DIR"/barcode03/fastq/cpu8 -o "$WORK_DIR"/barcode03/cpu8/search_output -s GTGATCTCGCTC -t 20
python "$CODE_DIR"/search_barcode/search_fastq_dyn_pro.py -i "$WORK_DIR"/barcode04/fastq/cpu8 -o "$WORK_DIR"/barcode04/cpu8/search_output -s GTGATCTCGCTC -t 20
python "$CODE_DIR"/search_barcode/search_fastq_dyn_pro.py -i "$WORK_DIR"/barcode05/fastq/cpu8 -o "$WORK_DIR"/barcode05/cpu8/search_output -s GTGATCTCGCTC -t 20
python "$CODE_DIR"/search_barcode/search_fastq_dyn_pro.py -i "$WORK_DIR"/barcode06/fastq/cpu8 -o "$WORK_DIR"/barcode06/cpu8/search_output -s GTGATCTCGCTC -t 20
python "$CODE_DIR"/search_barcode/search_fastq_dyn_pro.py -i "$WORK_DIR"/barcode07/fastq/cpu8 -o "$WORK_DIR"/barcode07/cpu8/search_output -s GTGATCTCGCTC -t 20
python "$CODE_DIR"/search_barcode/search_fastq_dyn_pro.py -i "$WORK_DIR"/barcode08/fastq/cpu8 -o "$WORK_DIR"/barcode08/cpu8/search_output -s GTGATCTCGCTC -t 20
python "$CODE_DIR"/search_barcode/search_fastq_dyn_pro.py -i "$WORK_DIR"/barcode09/fastq/cpu8 -o "$WORK_DIR"/barcode09/cpu8/search_output -s GTGATCTCGCTC -t 20
python "$CODE_DIR"/search_barcode/search_fastq_dyn_pro.py -i "$WORK_DIR"/barcode10/fastq/cpu8 -o "$WORK_DIR"/barcode10/cpu8/search_output -s GTGATCTCGCTC -t 20




python "$CODE_DIR"/search_barcode/search_barcode_dyn_pro.py -i "$WORK_DIR"/barcode01/cpu8/search_output -e "$CODE_DIR"/search_barcode/barcode_seq.xlsx -o "$WORK_DIR"/barcode01/cpu8/barcode_match1_output -s1 CCGCTGGTGGCT -s2 CACTCGTCGTGC  
python "$CODE_DIR"/search_barcode/search_barcode_dyn_pro.py -i "$WORK_DIR"/barcode01/cpu8/search_output -e "$CODE_DIR"/search_barcode/barcode_seq.xlsx -o "$WORK_DIR"/barcode01/cpu8/barcode_match2_output -s1 CCGCTGGTGGCT -s2 CACTCGTCGTGC  -r

python "$CODE_DIR"/search_barcode/search_barcode_dyn_pro.py -i "$WORK_DIR"/barcode02/cpu8/search_output -e "$CODE_DIR"/search_barcode/barcode_seq.xlsx -o "$WORK_DIR"/barcode02/cpu8/barcode_match1_output -s1 CCGCTGGTGGCT -s2 CACTCGTCGTGC  
python "$CODE_DIR"/search_barcode/search_barcode_dyn_pro.py -i "$WORK_DIR"/barcode02/cpu8/search_output -e "$CODE_DIR"/search_barcode/barcode_seq.xlsx -o "$WORK_DIR"/barcode02/cpu8/barcode_match2_output -s1 CCGCTGGTGGCT -s2 CACTCGTCGTGC  -r

python "$CODE_DIR"/search_barcode/search_barcode_dyn_pro.py -i "$WORK_DIR"/barcode03/cpu8/search_output -e "$CODE_DIR"/search_barcode/barcode_seq.xlsx -o "$WORK_DIR"/barcode03/cpu8/barcode_match1_output -s1 CCGCTGGTGGCT -s2 CACTCGTCGTGC  
python "$CODE_DIR"/search_barcode/search_barcode_dyn_pro.py -i "$WORK_DIR"/barcode03/cpu8/search_output -e "$CODE_DIR"/search_barcode/barcode_seq.xlsx -o "$WORK_DIR"/barcode03/cpu8/barcode_match2_output -s1 CCGCTGGTGGCT -s2 CACTCGTCGTGC  -r

python "$CODE_DIR"/search_barcode/search_barcode_dyn_pro.py -i "$WORK_DIR"/barcode04/cpu8/search_output -e "$CODE_DIR"/search_barcode/barcode_seq.xlsx -o "$WORK_DIR"/barcode04/cpu8/barcode_match1_output -s1 CCGCTGGTGGCT -s2 CACTCGTCGTGC  
python "$CODE_DIR"/search_barcode/search_barcode_dyn_pro.py -i "$WORK_DIR"/barcode04/cpu8/search_output -e "$CODE_DIR"/search_barcode/barcode_seq.xlsx -o "$WORK_DIR"/barcode04/cpu8/barcode_match2_output -s1 CCGCTGGTGGCT -s2 CACTCGTCGTGC  -r

python "$CODE_DIR"/search_barcode/search_barcode_dyn_pro.py -i "$WORK_DIR"/barcode05/cpu8/search_output -e "$CODE_DIR"/search_barcode/barcode_seq.xlsx -o "$WORK_DIR"/barcode05/cpu8/barcode_match1_output -s1 CCGCTGGTGGCT -s2 CACTCGTCGTGC  
python "$CODE_DIR"/search_barcode/search_barcode_dyn_pro.py -i "$WORK_DIR"/barcode05/cpu8/search_output -e "$CODE_DIR"/search_barcode/barcode_seq.xlsx -o "$WORK_DIR"/barcode05/cpu8/barcode_match2_output -s1 CCGCTGGTGGCT -s2 CACTCGTCGTGC  -r

python "$CODE_DIR"/search_barcode/search_barcode_dyn_pro.py -i "$WORK_DIR"/barcode06/cpu8/search_output -e "$CODE_DIR"/search_barcode/barcode_seq.xlsx -o "$WORK_DIR"/barcode06/cpu8/barcode_match1_output -s1 CCGCTGGTGGCT -s2 CACTCGTCGTGC  
python "$CODE_DIR"/search_barcode/search_barcode_dyn_pro.py -i "$WORK_DIR"/barcode06/cpu8/search_output -e "$CODE_DIR"/search_barcode/barcode_seq.xlsx -o "$WORK_DIR"/barcode06/cpu8/barcode_match2_output -s1 CCGCTGGTGGCT -s2 CACTCGTCGTGC  -r

python "$CODE_DIR"/search_barcode/search_barcode_dyn_pro.py -i "$WORK_DIR"/barcode07/cpu8/search_output -e "$CODE_DIR"/search_barcode/barcode_seq.xlsx -o "$WORK_DIR"/barcode07/cpu8/barcode_match1_output -s1 CCGCTGGTGGCT -s2 CACTCGTCGTGC  
python "$CODE_DIR"/search_barcode/search_barcode_dyn_pro.py -i "$WORK_DIR"/barcode07/cpu8/search_output -e "$CODE_DIR"/search_barcode/barcode_seq.xlsx -o "$WORK_DIR"/barcode07/cpu8/barcode_match2_output -s1 CCGCTGGTGGCT -s2 CACTCGTCGTGC  -r

python "$CODE_DIR"/search_barcode/search_barcode_dyn_pro.py -i "$WORK_DIR"/barcode08/cpu8/search_output -e "$CODE_DIR"/search_barcode/barcode_seq.xlsx -o "$WORK_DIR"/barcode08/cpu8/barcode_match1_output -s1 CCGCTGGTGGCT -s2 CACTCGTCGTGC  
python "$CODE_DIR"/search_barcode/search_barcode_dyn_pro.py -i "$WORK_DIR"/barcode08/cpu8/search_output -e "$CODE_DIR"/search_barcode/barcode_seq.xlsx -o "$WORK_DIR"/barcode08/cpu8/barcode_match2_output -s1 CCGCTGGTGGCT -s2 CACTCGTCGTGC  -r

python "$CODE_DIR"/search_barcode/search_barcode_dyn_pro.py -i "$WORK_DIR"/barcode09/cpu8/search_output -e "$CODE_DIR"/search_barcode/barcode_seq.xlsx -o "$WORK_DIR"/barcode09/cpu8/barcode_match1_output -s1 CCGCTGGTGGCT -s2 CACTCGTCGTGC  
python "$CODE_DIR"/search_barcode/search_barcode_dyn_pro.py -i "$WORK_DIR"/barcode09/cpu8/search_output -e "$CODE_DIR"/search_barcode/barcode_seq.xlsx -o "$WORK_DIR"/barcode09/cpu8/barcode_match2_output -s1 CCGCTGGTGGCT -s2 CACTCGTCGTGC  -r

python "$CODE_DIR"/search_barcode/search_barcode_dyn_pro.py -i "$WORK_DIR"/barcode10/cpu8/search_output -e "$CODE_DIR"/search_barcode/barcode_seq.xlsx -o "$WORK_DIR"/barcode10/cpu8/barcode_match1_output -s1 CCGCTGGTGGCT -s2 CACTCGTCGTGC  
python "$CODE_DIR"/search_barcode/search_barcode_dyn_pro.py -i "$WORK_DIR"/barcode10/cpu8/search_output -e "$CODE_DIR"/search_barcode/barcode_seq.xlsx -o "$WORK_DIR"/barcode10/cpu8/barcode_match2_output -s1 CCGCTGGTGGCT -s2 CACTCGTCGTGC  -r
