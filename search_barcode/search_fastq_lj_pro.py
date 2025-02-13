#!/bin/bash 
# -*- coding: utf-8 -*-
import os
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Align import PairwiseAligner
import multiprocessing as mp


# 设置参数解析器
parser = argparse.ArgumentParser(description='Find a specific sequence in FASTQ files')
parser.add_argument('-i', '--fastq_folder', required=True, type=str, help='Path to the folder containing FASTQ files')
parser.add_argument('-o', '--output_folder', required=True, help='Output fastq file for matched sequences')
parser.add_argument('-s', '--specific_seq', required=True, type=str, help='The specific sequence to be found in the FASTQ files')
parser.add_argument('-t', '--threshold', type=int, default=20, help='The threshold score for matching sequences (default=20)')
args = parser.parse_args()

if not os.path.exists(args.output_folder):
    os.makedirs(args.output_folder)

# 定义函数：反向互补序列
def reverse_complement(seq):
    return str(Seq(seq).reverse_complement())


# 定义函数：Smith-Waterman算法匹配
def smith_waterman(seq1, seq2):
    # 初始化比对器
    aligner = PairwiseAligner()
    aligner.mode = 'local'  # 设置为局部比对
    aligner.match_score = 2  # 匹配得分
    aligner.mismatch_score = -1  # 不匹配得分
    aligner.open_gap_score = -1  # 打开缺口罚分
    aligner.extend_gap_score = -1  # 扩展缺口罚分

    # 执行比对
    alignments = aligner.align(seq1, seq2)

    # 如果没有比对结果，返回默认值
    if not alignments:
        return [None, None, 0, 0, None, None]

    # 按得分排序比对结果
    sorted_alignments = sorted(alignments, key=lambda x: -x.score)

    # 获取最佳比对结果
    top_alignment = sorted_alignments[0]
    top_score = top_alignment.score
    top_aln_start = top_alignment.aligned[0][0][0]  # 比对起始位置
    top_aln_end = top_alignment.aligned[0][-1][-1]  # 比对结束位置

    # 获取次优比对得分
    try:
        second_top_score = sorted_alignments[1].score
    except IndexError:
        second_top_score = 0

    # 获取比对序列
    seq1_match = str(top_alignment.target)
    seq2_match = str(top_alignment.query)

    return [seq1_match, seq2_match, int(top_score), int(second_top_score), int(top_aln_start), int(top_aln_end)]

def process_fastq(file_path, specific_seq, threshold, process_id, output_folder):
    # 创建每个进程的临时文件夹
    temp_folder = os.path.join(output_folder, f"process_{process_id}")
    os.makedirs(temp_folder, exist_ok=True)
    output_path_match1 = os.path.join(temp_folder, "output_match1.fastq")
    output_path_match2 = os.path.join(temp_folder, "output_match2.fastq")

    total_seq = 0
    matched_seq = 0
    re_matched_seq = 0
    complementary_matched_seq = 0
    multiple_matched_seq = 0
    fail_matched_seq = 0
    results = []

    for record in SeqIO.parse(file_path, "fastq"):
        seq = str(record.seq)

        # 匹配正向序列
        top_aln = smith_waterman(seq, specific_seq)
        # 匹配反向互补序列
        rc_top_aln = smith_waterman(seq, reverse_complement(specific_seq))

        # 判断匹配是否成功
        matched = 0 #1代表正向匹配，2代表反向匹配，3代表存在互补barcode，4代表多正或多反barcode，5代表没有barcode
        if top_aln[2] >= threshold and rc_top_aln[2] >= threshold:
            complementary_matched_seq += 1
        elif top_aln[3] >= threshold or rc_top_aln[3] >= threshold:
            multiple_matched_seq += 1
        elif top_aln[2] < threshold and rc_top_aln[2] < threshold:
            fail_matched_seq += 1
        elif top_aln[2] >= threshold:
            matched_seq += 1
            matched = 1
            sub_seq = top_aln[0][top_aln[4]-29:top_aln[4]]
            results.append((record, matched, sub_seq))

        elif rc_top_aln[2] >= threshold:
            re_matched_seq += 1
            matched = 2
            sub_seq = rc_top_aln[0][rc_top_aln[5]:rc_top_aln[5]+29]
            results.append((record, matched, sub_seq))
        
        total_seq += 1
    
    with open(output_path_match1, "a") as out_file1, open(output_path_match2, "a") as out_file2:
        for record, match, sub_seq in results:
            if match == 1:
                write_output(record, sub_seq, out_file1)
            elif match == 2:
                write_output(record, sub_seq, out_file2)

    return [total_seq, matched_seq, re_matched_seq,complementary_matched_seq, multiple_matched_seq, fail_matched_seq, temp_folder]

def write_output(record, sub_seq, out_file):
    seq = str(record.seq)
    out_file.write('@' + record.id + '_'+str(sub_seq)+'\n')
    out_file.write(seq + '\n')
    out_file.write('+' + '\n')
    quality_string = ''.join([chr(score + 33) for score in record.letter_annotations['phred_quality']])
    out_file.write(quality_string + '\n')

#def merge_and_clean(temp_folders, final_output_paths):
#    for final_output_path in final_output_paths:
#        with open(final_output_path, "w") as final_file:
#            for folder in temp_folders:
#                temp_file_path = os.path.join(folder, os.path.basename(final_output_path))
#                if os.path.exists(temp_file_path):
#                    with open(temp_file_path, "r") as temp_file:
#                        shutil.copyfileobj(temp_file, final_file)
#    
#    # 删除临时文件夹，这一步在所有必要文件合并后执行
#    for folder in temp_folders:
#        if os.path.exists(folder):
#            shutil.rmtree(folder)
      

def main():
    fastq_files = [os.path.join(args.fastq_folder, file_name) for file_name in os.listdir(args.fastq_folder) if file_name.endswith('.fastq')]

    pool = mp.Pool(mp.cpu_count())

    # 使用 starmap 时，传递所有必要的参数
    args_for_starmap = [(file_path, args.specific_seq, args.threshold, process_id, args.output_folder) for process_id, file_path in enumerate(fastq_files)]
    results = pool.starmap(process_fastq, args_for_starmap)

    # 合并临时文件并清理
    #temp_folders = [os.path.join(args.output_folder, f"process_{pid}") for pid in range(len(fastq_files))]
    #final_output_paths = [os.path.join(args.output_folder, "output_match1.fastq"), os.path.join(args.output_folder, "output_match2.fastq")]
    #merge_and_clean(temp_folders, final_output_paths)

    # 输出统计结果
    total_seq = sum(res[0] for res in results)
    matched_seq = sum(res[1] for res in results)
    re_matched_seq = sum(res[2] for res in results)
    complementary_matched_seq = sum(res[3] for res in results)
    multiple_matched_seq = sum(res[4] for res in results)
    fail_matched_seq = sum(res[5] for res in results)
    parent_folder = os.path.dirname(args.output_folder)
    output_file_path = os.path.join(parent_folder, "statistical.txt")
    with open(output_file_path, "w") as f:
        f.write('序列总数：' + str(total_seq) + '\n')
        f.write('单次正向匹配序列百分比：' + str(matched_seq / total_seq) + '\n')
        f.write('单次反向匹配序列百分比：' + str(re_matched_seq / total_seq) + '\n')
        f.write('互补匹配序列百分比：' + str(complementary_matched_seq / total_seq) + '\n')
        f.write('多次匹配序列百分比：' + str(multiple_matched_seq / total_seq) + '\n')
        f.write('匹配失败序列百分比：' + str(fail_matched_seq / total_seq) + '\n')

if __name__ == "__main__":
    main()


