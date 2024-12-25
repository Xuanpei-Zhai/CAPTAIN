import os
import shutil
import argparse
from Bio import SeqIO
from collections import defaultdict
import matplotlib.pyplot as plt
def count_sequences_in_fastq_files(input_directory, output_file_path, threshold):
    sequence_count_dict = defaultdict(int)
    total_files = 0
    moved_files = 0
    total_sequences_before_moving = 0
    total_sequences_after_moving = 0
    sequence_counts = []

    # Create a new directory to store moved files
    new_directory = os.path.join(output_file_path, 'bbmap_output_many_seqs')
    if not os.path.exists(new_directory):
        os.makedirs(new_directory)

    # Open a log file to record errors
    with open(os.path.join(output_file_path, "error_log.txt"), 'w') as error_log:
        for filename in os.listdir(input_directory):
            if filename.endswith(".fastq"):
                total_files += 1
                file_path = os.path.join(input_directory, filename)
                sequence_count = 0
                try:
                    for _ in SeqIO.parse(file_path, "fastq"):
                        sequence_count += 1
                except ValueError as e:
                    # Log the error and continue with the next file
                    error_msg = f"Error processing {filename}: {e}\n"
                    print(error_msg)  # Print to console or consider logging to a file
                    error_log.write(error_msg)
                    continue  # Skip the rest of the current iteration and proceed with the next file

                total_sequences_before_moving += sequence_count

                if sequence_count >= threshold:
                    shutil.copy(file_path, new_directory)  # Move the file instead of deleting
                    moved_files += 1
                else:
                    total_sequences_after_moving += sequence_count
                sequence_count_dict[sequence_count] += 1
                sequence_counts.append(sequence_count)

        sequence_counts.sort(reverse=True)

        with open(os.path.join(output_file_path, "bbmap_statistical.txt"), 'w') as output_file:
            output_file.write(f"Total files: {total_files}\n")
            output_file.write(f"Moved files: {moved_files}\n")
            output_file.write(f"Total sequences before moving: {total_sequences_before_moving}\n")
            output_file.write(f"Total sequences after moving: {total_sequences_after_moving}\n")
            output_file.write("Number of files with specific sequence counts:\n")
            for sequence_count, file_count in sorted(sequence_count_dict.items()):
                output_file.write(f"{sequence_count} sequences: {file_count} files\n")
        # Create histogram and curve fitting
        x_values = range(1, len(sequence_counts) + 1)
        plt.plot(x_values,sequence_counts, marker='o', linestyle='-')
        plt.yscale("log")
        plt.xscale("log")
        plt.xlabel("barcodes(log scale)")
        plt.ylabel("Number of sequences (log scale)")
        plt.title("Fastq files sequence count distribution (sequence count > 1)")

        # Set x and y axis limits
        #plt.xlim(1, 60000)
        #plt.ylim(1, 1000)


        plt.savefig(args.output_file_path + "/histogram_and_curve_gt_1_log_scale.png", format='png')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process and filter fastq files")
    parser.add_argument('-i', '--input_directory', required=True, type=str, help="Input directory path")
    parser.add_argument('-o', '--output_file_path', required=True, type=str, help="Output file path")
    parser.add_argument('-s', '--threshold', required=True, type=int, help="Sequence count threshold")

    args = parser.parse_args()
    count_sequences_in_fastq_files(args.input_directory, args.output_file_path, args.threshold)
