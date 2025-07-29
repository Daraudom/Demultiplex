#!/usr/bin/env python
import numpy as np
import concurrent.futures
import argparse
import gzip
import matplotlib.pyplot as plt
import sys


# define function to parse cmd-line arguments
def get_args():
    parser = argparse.ArgumentParser(description="Generates a quality score histogram from\
                                      a forward read fq file, forward barcode fq file,\
                                     reverse read fq file and reverse barcode file")
    parser.add_argument('-r1', '--read1', help="Give a forward read fq file.")
    parser.add_argument('-r2', '--read2', help="Give a forward index fq file.")
    parser.add_argument('-r3', '--read3', help="Give a reverse read fq file.")
    parser.add_argument('-r4', '--read4', help="Give a reverse index fq file.")
    return parser.parse_args()

# Establish Global Variables
args = get_args()
read1 = args.read1
read2 = args.read2
read3 = args.read3
read4 = args.read4
records_num = 1452986940/4
read_length = 101

# Outputs the quality score distribution of the reads
def process_read(read_file):
    """
    Input: fastq file
    """
    # initialize an dictionary-> {base_position: [list of qual scores]}
    #bp_scores = {i : [] for i in range(101)}
    try:
        bp_scores = {}

        with gzip.open(read_file, 'rb') as fin:
            while True:
                header = fin.readline().strip()
                # break loop
                if header == b'': break
                seq = fin.readline().strip()
                plus = fin.readline().strip()
                qual = fin.readline().strip()

                if len(bp_scores) == 0:
                    bp_scores = {i : 0 for i in range(len(seq))}

                # Iterate through the quality strings
                for i in range(len(qual)):
                    bp_scores[i]+=(int(qual[i])-33)
        
        # Turn the list value into an np array
        for bp, qual in bp_scores.items():
            bp_scores[bp] = qual/records_num # type: ignore

    except Exception as e:
        print(f"Error processing {read_file}: {e}", file=sys.stderr)
        return {}

    return bp_scores

# Plot the histogram
def plot_distribution(read_dict, label=1):
    y_axis = list(read_dict.values())
    x_axis = np.arange(0, len(read_dict))

    plt.bar(x_axis, y_axis, width=0.5)
    plt.xlabel('Base Position')
    plt.ylabel('Mean Quality Score')
    plt.title(f'Mean Quality Score Distribution of Read {label}')
    plt.grid(visible=True, which='both', alpha=0.5)
    plt.savefig(f'qual_score_distribution_read{label}.png')
    plt.close()

reads = [read1, read2, read3, read4]

with concurrent.futures.ProcessPoolExecutor() as executor:
    results = list(executor.map(process_read, reads))
    print(results)
    
    for i, bp_dt in enumerate(results, start=1):
        print(bp_dt)
        plot_distribution(bp_dt, label=i)
        



