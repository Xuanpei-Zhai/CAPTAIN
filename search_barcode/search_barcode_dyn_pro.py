#!/bin/bash 
import argparse
import pandas as pd
import os
import multiprocessing as mp
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Align import PairwiseAligner
import glob


def reverse_complement(sequence):
    return str(Seq(sequence).reverse_complement())


def fuzzy_match_global(query_seq, white_list):
    aligner = PairwiseAligner()
    aligner.mode = 'global'
    aligner.match_score = 1
    aligner.mismatch_score = -1
    aligner.open_gap_score = -1
    aligner.extend_gap_score = -1

    for subseq in white_list:
        alignments = aligner.align(subseq, query_seq)
        if alignments and alignments[0].score >= len(subseq) - 2:
            query_seq = subseq
            return True, query_seq
    return False, query_seq

def fuzzy_match_local(query_seq, white_list):
    aligner = PairwiseAligner()
    aligner.mode = 'local'
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -1
    aligner.extend_gap_score = -1

    for subseq in white_list:
        alignments = aligner.align(subseq, query_seq)
        if alignments and alignments[0].score >= 2 * len(subseq) - 7:
            match_start = alignments[0].aligned[0][0][0]
            match_end = alignments[0].aligned[0][-1][-1]
            return True, query_seq[match_start:match_end], match_start, match_end
    return False, query_seq, None, None


def process_file(fastq_file, col1, col2, col3, output_folder, fixed_seq1, fixed_seq2, reverse, process_id):

    # 创建每个进程的子文件夹
    process_output_folder = os.path.join(output_folder, f"process_{process_id}")
    os.makedirs(process_output_folder, exist_ok=True)
    output_path = os.path.join(process_output_folder, f"{os.path.basename(fastq_file)}.fastq")
    
    counts = {'total_seqs': 0, 'link_right': 0, 'barcode1_count': 0, 'barcode2_count': 0, 'barcode3_count': 0, 'triple_barcodes_count': 0}
    results = []

    with open(output_path, "w") as output_file: 
        for record in SeqIO.parse(fastq_file, "fastq"):
            counts['total_seqs'] += 1

            barcode_sequence = record.description.split("_")[-1]
            if reverse:
                sequence = reverse_complement(barcode_sequence)
            else:
                sequence = barcode_sequence

            fuzzy_match_seq1 = fuzzy_match_local(sequence, [fixed_seq1])
            fuzzy_match_seq2 = fuzzy_match_local(sequence, [fixed_seq2])
            if fuzzy_match_seq1[0] and fuzzy_match_seq2[0]:
                index1_start = fuzzy_match_seq1[2]
                index1_end = index1_start + len(fuzzy_match_seq1[1])
                index2_start = fuzzy_match_seq2[2]
                index2_end = index2_start + len(fuzzy_match_seq2[1])

                if index1_end <= index2_start:
                    counts['link_right'] += 1

                    barcode1 = sequence[:index1_start]
                    barcode2 = sequence[index1_end:index2_start]
                    barcode3 = sequence[index2_end:]


                    barcode1_flag = 0
                    barcode2_flag = 0
                    barcode3_flag = 0

                    fuzzy_match_barcode1 = fuzzy_match_global(barcode1, col3)
                    fuzzy_match_barcode2 = fuzzy_match_global(barcode2, col2)
                    fuzzy_match_barcode3 = fuzzy_match_global(barcode3, col1)

                    if fuzzy_match_barcode1[0]:
                    # if barcode1 in col1 or reverse_complement(barcode1) in col3:
                        barcode1_flag = 1
                        counts['barcode1_count'] += 1
                        barcode1 = fuzzy_match_barcode1[1]

                    if fuzzy_match_barcode2[0]:
                    # if barcode2 in col2 or reverse_complement(barcode2) in col2:
                        barcode2_flag = 1
                        counts['barcode2_count'] += 1
                        barcode2 = fuzzy_match_barcode2[1]
                    if fuzzy_match_barcode3[0]:
                    # if barcode3 in col3 or reverse_complement(barcode3) in col1:
                        barcode3_flag = 1
                        counts['barcode3_count'] += 1
                        barcode3 = fuzzy_match_barcode3[1]
                    if barcode1_flag and barcode2_flag and barcode3_flag:
                        counts['triple_barcodes_count'] += 1
                        barcode_name = barcode1 +fixed_seq1+ barcode2 +fixed_seq2+barcode3
                        output_file.write('@' + record.id.split('_')[0] + '_'+str(barcode_name)+'\n')
                        output_file.write(str(record.seq) + '\n')
                        output_file.write('+' + '\n')
                        quality_string = ''.join([chr(score + 33) for score in record.letter_annotations['phred_quality']])
                        output_file.write(quality_string + '\n')

        return counts

def count_sequences(directory, excel_file, output_folder, fixed_seq1, fixed_seq2, reverse):

    data = pd.read_excel(excel_file)
    col1, col2, col3 = data.iloc[:, 0].tolist(), data.iloc[:, 1].tolist(), data.iloc[:, 2].tolist()


    if reverse:
        fastq_files = glob.glob(os.path.join(directory, '**/output_match2.fastq'), recursive=True)
    else:
        fastq_files = glob.glob(os.path.join(directory, '**/output_match1.fastq'), recursive=True)


    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    pool = mp.Pool(mp.cpu_count())

    args_for_starmap = [(fastq_file, col1, col2, col3, output_folder, fixed_seq1, fixed_seq2, reverse, process_id)
                        for process_id, fastq_file in enumerate(fastq_files)]

    results = pool.starmap(process_file, args_for_starmap)
    pool.close()
    pool.join()


    # 汇总统计信息
    summary = {'total_seqs': 0, 'link_right': 0, 'barcode1_count': 0, 'barcode2_count': 0, 'barcode3_count': 0, 'triple_barcodes_count': 0}
    for counts in results:
        for key in summary:
            summary[key] += counts[key]

    with open(os.path.join(output_folder, "summary.txt"), "w") as summary_file:
        for key, value in summary.items():
            summary_file.write(f"{key}: {value}\n")



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process fasta and excel files")
    parser.add_argument('-i',"--directory",required=True,type=str, help="Directory path")
    parser.add_argument('-e',"--excel_file",required=True, type=str,help="Excel file path")
    parser.add_argument('-o', "--output_folder", required=True, type=str, help="Output folder path")
    parser.add_argument('-s1',"--fixed_seq1",required=True,type=str, help="Fixed sequence 1")
    parser.add_argument('-s2',"--fixed_seq2",required=True,type=str, help="Fixed sequence 2")
    parser.add_argument('-r', "--reverse", required=False, default=False, action='store_true', help="Flag to indicate if a reverse complement sequence is needed")

    args = parser.parse_args()

    count_sequences(args.directory, args.excel_file,args.output_folder, args.fixed_seq1, args.fixed_seq2, args.reverse)
