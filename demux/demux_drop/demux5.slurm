#!/bin/bash
#SBATCH -J step3_5_demux
#SBATCH -p dm_pub_cpu
#SBATCH -o /dm_data/zhaixp/step3_5_demux_out
#SBATCH -e /dm_data/zhaixp/step3_5_demux_error
#SBATCH -N 1
#SBATCH -n 32
#SBATCH --time=3-00:00:00
#SBATCH --exclude=dm-c02n03

source ~/CAPTAIN/workflow/config.sh
source ~/anaconda3/etc/profile.d/conda.sh
conda activate captain

python "$CODE_DIR"/search_barcode/search_fastq_lj_pro.py -i "$WORK_DIR"/barcode01/fastq/cpu5 -o "$WORK_DIR"/barcode01/cpu5/search_output -s GTCTCGTGGGCTCGG -t 26
python "$CODE_DIR"/search_barcode/search_fastq_lj_pro.py -i "$WORK_DIR"/barcode02/fastq/cpu5 -o "$WORK_DIR"/barcode02/cpu5/search_output -s GTCTCGTGGGCTCGG -t 26
python "$CODE_DIR"/search_barcode/search_fastq_lj_pro.py -i "$WORK_DIR"/barcode03/fastq/cpu5 -o "$WORK_DIR"/barcode03/cpu5/search_output -s GTCTCGTGGGCTCGG -t 26
python "$CODE_DIR"/search_barcode/search_fastq_lj_pro.py -i "$WORK_DIR"/barcode04/fastq/cpu5 -o "$WORK_DIR"/barcode04/cpu5/search_output -s GTCTCGTGGGCTCGG -t 26
python "$CODE_DIR"/search_barcode/search_fastq_lj_pro.py -i "$WORK_DIR"/barcode05/fastq/cpu5 -o "$WORK_DIR"/barcode05/cpu5/search_output -s GTCTCGTGGGCTCGG -t 26
python "$CODE_DIR"/search_barcode/search_fastq_lj_pro.py -i "$WORK_DIR"/barcode06/fastq/cpu5 -o "$WORK_DIR"/barcode06/cpu5/search_output -s GTCTCGTGGGCTCGG -t 26
python "$CODE_DIR"/search_barcode/search_fastq_lj_pro.py -i "$WORK_DIR"/barcode07/fastq/cpu5 -o "$WORK_DIR"/barcode07/cpu5/search_output -s GTCTCGTGGGCTCGG -t 26
python "$CODE_DIR"/search_barcode/search_fastq_lj_pro.py -i "$WORK_DIR"/barcode08/fastq/cpu5 -o "$WORK_DIR"/barcode08/cpu5/search_output -s GTCTCGTGGGCTCGG -t 26
python "$CODE_DIR"/search_barcode/search_fastq_lj_pro.py -i "$WORK_DIR"/barcode09/fastq/cpu5 -o "$WORK_DIR"/barcode09/cpu5/search_output -s GTCTCGTGGGCTCGG -t 26
python "$CODE_DIR"/search_barcode/search_fastq_lj_pro.py -i "$WORK_DIR"/barcode10/fastq/cpu5 -o "$WORK_DIR"/barcode10/cpu5/search_output -s GTCTCGTGGGCTCGG -t 26




python "$CODE_DIR"/search_barcode/search_barcode_lj_pro.py -i "$WORK_DIR"/barcode01/cpu5/search_output -e "$CODE_DIR"/search_barcode/barcode_seq.xlsx -o "$WORK_DIR"/barcode01/cpu5/barcode_match1_output -s1 ACAG -s2 GTCA 
python "$CODE_DIR"/search_barcode/search_barcode_lj_pro.py -i "$WORK_DIR"/barcode01/cpu5/search_output -e "$CODE_DIR"/search_barcode/barcode_seq.xlsx -o "$WORK_DIR"/barcode01/cpu5/barcode_match2_output -s1 ACAG -s2 GTCA -r

python "$CODE_DIR"/search_barcode/search_barcode_lj_pro.py -i "$WORK_DIR"/barcode02/cpu5/search_output -e "$CODE_DIR"/search_barcode/barcode_seq.xlsx -o "$WORK_DIR"/barcode02/cpu5/barcode_match1_output -s1 ACAG -s2 GTCA 
python "$CODE_DIR"/search_barcode/search_barcode_lj_pro.py -i "$WORK_DIR"/barcode02/cpu5/search_output -e "$CODE_DIR"/search_barcode/barcode_seq.xlsx -o "$WORK_DIR"/barcode02/cpu5/barcode_match2_output -s1 ACAG -s2 GTCA -r

python "$CODE_DIR"/search_barcode/search_barcode_lj_pro.py -i "$WORK_DIR"/barcode03/cpu5/search_output -e "$CODE_DIR"/search_barcode/barcode_seq.xlsx -o "$WORK_DIR"/barcode03/cpu5/barcode_match1_output -s1 ACAG -s2 GTCA 
python "$CODE_DIR"/search_barcode/search_barcode_lj_pro.py -i "$WORK_DIR"/barcode03/cpu5/search_output -e "$CODE_DIR"/search_barcode/barcode_seq.xlsx -o "$WORK_DIR"/barcode03/cpu5/barcode_match2_output -s1 ACAG -s2 GTCA -r

python "$CODE_DIR"/search_barcode/search_barcode_lj_pro.py -i "$WORK_DIR"/barcode04/cpu5/search_output -e "$CODE_DIR"/search_barcode/barcode_seq.xlsx -o "$WORK_DIR"/barcode04/cpu5/barcode_match1_output -s1 ACAG -s2 GTCA 
python "$CODE_DIR"/search_barcode/search_barcode_lj_pro.py -i "$WORK_DIR"/barcode04/cpu5/search_output -e "$CODE_DIR"/search_barcode/barcode_seq.xlsx -o "$WORK_DIR"/barcode04/cpu5/barcode_match2_output -s1 ACAG -s2 GTCA -r

python "$CODE_DIR"/search_barcode/search_barcode_lj_pro.py -i "$WORK_DIR"/barcode05/cpu5/search_output -e "$CODE_DIR"/search_barcode/barcode_seq.xlsx -o "$WORK_DIR"/barcode05/cpu5/barcode_match1_output -s1 ACAG -s2 GTCA 
python "$CODE_DIR"/search_barcode/search_barcode_lj_pro.py -i "$WORK_DIR"/barcode05/cpu5/search_output -e "$CODE_DIR"/search_barcode/barcode_seq.xlsx -o "$WORK_DIR"/barcode05/cpu5/barcode_match2_output -s1 ACAG -s2 GTCA -r

python "$CODE_DIR"/search_barcode/search_barcode_lj_pro.py -i "$WORK_DIR"/barcode06/cpu5/search_output -e "$CODE_DIR"/search_barcode/barcode_seq.xlsx -o "$WORK_DIR"/barcode06/cpu5/barcode_match1_output -s1 ACAG -s2 GTCA 
python "$CODE_DIR"/search_barcode/search_barcode_lj_pro.py -i "$WORK_DIR"/barcode06/cpu5/search_output -e "$CODE_DIR"/search_barcode/barcode_seq.xlsx -o "$WORK_DIR"/barcode06/cpu5/barcode_match2_output -s1 ACAG -s2 GTCA -r

python "$CODE_DIR"/search_barcode/search_barcode_lj_pro.py -i "$WORK_DIR"/barcode07/cpu5/search_output -e "$CODE_DIR"/search_barcode/barcode_seq.xlsx -o "$WORK_DIR"/barcode07/cpu5/barcode_match1_output -s1 ACAG -s2 GTCA 
python "$CODE_DIR"/search_barcode/search_barcode_lj_pro.py -i "$WORK_DIR"/barcode07/cpu5/search_output -e "$CODE_DIR"/search_barcode/barcode_seq.xlsx -o "$WORK_DIR"/barcode07/cpu5/barcode_match2_output -s1 ACAG -s2 GTCA -r

python "$CODE_DIR"/search_barcode/search_barcode_lj_pro.py -i "$WORK_DIR"/barcode08/cpu5/search_output -e "$CODE_DIR"/search_barcode/barcode_seq.xlsx -o "$WORK_DIR"/barcode08/cpu5/barcode_match1_output -s1 ACAG -s2 GTCA 
python "$CODE_DIR"/search_barcode/search_barcode_lj_pro.py -i "$WORK_DIR"/barcode08/cpu5/search_output -e "$CODE_DIR"/search_barcode/barcode_seq.xlsx -o "$WORK_DIR"/barcode08/cpu5/barcode_match2_output -s1 ACAG -s2 GTCA -r

python "$CODE_DIR"/search_barcode/search_barcode_lj_pro.py -i "$WORK_DIR"/barcode09/cpu5/search_output -e "$CODE_DIR"/search_barcode/barcode_seq.xlsx -o "$WORK_DIR"/barcode09/cpu5/barcode_match1_output -s1 ACAG -s2 GTCA 
python "$CODE_DIR"/search_barcode/search_barcode_lj_pro.py -i "$WORK_DIR"/barcode09/cpu5/search_output -e "$CODE_DIR"/search_barcode/barcode_seq.xlsx -o "$WORK_DIR"/barcode09/cpu5/barcode_match2_output -s1 ACAG -s2 GTCA -r

python "$CODE_DIR"/search_barcode/search_barcode_lj_pro.py -i "$WORK_DIR"/barcode10/cpu5/search_output -e "$CODE_DIR"/search_barcode/barcode_seq.xlsx -o "$WORK_DIR"/barcode10/cpu5/barcode_match1_output -s1 ACAG -s2 GTCA 
python "$CODE_DIR"/search_barcode/search_barcode_lj_pro.py -i "$WORK_DIR"/barcode10/cpu5/search_output -e "$CODE_DIR"/search_barcode/barcode_seq.xlsx -o "$WORK_DIR"/barcode10/cpu5/barcode_match2_output -s1 ACAG -s2 GTCA -r
