#!/usr/bin/env python

import argparse
import bioinfo as bf
import gzip
from tqdm import tqdm

# define a function to parse the cmd-line arguments

def get_args():
    parser = argparse.ArgumentParser(description="Given a pair of paired end reads (read1, index1\
                                     read2, index2), demultiplex the fq files.")
    read_group = parser.add_argument_group(title="Fastq Files StdIn", description='Input paired end biological reads\
                              and their index files.')
    read_group.add_argument('-r1', '--read1', help='Input the biological read 1.', required=True)
    read_group.add_argument('-r2', '--index1', help='Input the index read 1.', required= True)
    read_group.add_argument('-r3', '--index2', help='Input the index read 2.', required = True)
    read_group.add_argument('-r4', '--read2', help='Input the biological read 2.', required = True)
    read_group.add_argument('-idx', '--index_file', help='Input the list of indexes', required = True)

    qual_group = parser.add_argument_group(title="Quality Threshold Control", \
                                           description='Input quality score threshold for reads\
                          and index reads.')
    #qual_group.add_argument('-qr', '--read_threshold', default=35, required=False, help='Default is 35.')
    qual_group.add_argument('-qi', '--index_threshold', default=30, required=False, help='Default is 30.')

    return parser.parse_args()

# Allocating the arguments as a global variable int he script
args = get_args()
read1 = args.read1
index1 = args.index1
read2 = args.read2
index2 = args.index2
index_file = args.index_file
#qual_read = args.read_threshold
qual_index = int(args.index_threshold)

# Helper functions
def create_index_files(idx_file =index_file): 
    """
    Take an index file and return a list of index and a dictionary of index file handlers.
    Returns 52 writeable files one for R1 and R2.
    """

    idx_list = []
    idx_files_dict = {}

    with open(idx_file, 'r') as fin:
        # Skip past the header
        fin.readline()

        for line in fin:
            idx_list.append(line.split()[4])
    
    for idx in idx_list:
        idx_files_dict[idx] = []
        idx_files_dict[idx].append(open(f'{idx}_R1.fq', 'w'))
        idx_files_dict[idx].append(open(f'{idx}_R2.fq', 'w'))
    
    idx_files_dict['unk'] = [open('Unk_R1.fq', 'w'), open('Unk_R2.fq', 'w')]
    idx_files_dict['hop'] = [open('Hop_R1.fq', 'w'), open('Hop_R2.fq', 'w')]

    return idx_list, idx_files_dict

def close_idx_files(idx_dict:dict):
    """
    Close all Index Files Open.
    """
    for idx in idx_dict:
        print(idx)
        idx_dict[idx][0].close()
        idx_dict[idx][1].close()

def reverse_complement(seq:str):
    """
    Return the reverse complement of the sequence.
    """
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    comp = ''

    for char in seq:
        comp += complement[char]
    return comp[::-1]

def search_for_unks(idx1:str, idx2:str, idx1_qual:str, idx2_qual:str, idx_list: list[str], th: float) -> bool:
    """
    Check if index1 or index2 contain 'N' OR
    isn't found in idx_list OR
    is below the set quality idx threshold.
    """

    # Checking if there are uknowns
    if 'N' in idx1 or 'N' in idx2:
        return True
    
    # Checking if it is a valid index
    if idx1 not in idx_list or idx2 not in idx_list:
        return True
    
    # Checking Quality Scores
    for i in range(8):

        idx1_score = bf.convert_phred(idx1_qual[i])
        idx2_score = bf.convert_phred(idx2_qual[i])

        if idx1_score < th or idx2_score < th:
            return True
    
    return False

def search_for_hops(idx1:str, idx2:str) -> bool:
    """
    Check if index1 or index2 are swapped
    """
    if idx1 != idx2:
        return True
    return False

def read_fq(fq_handle):
    """
    Given a fq file handle, read the record.
    """
    read_header = fq_handle.readline().strip()
    read_seq = fq_handle.readline().strip()
    read_plus = fq_handle.readline().strip()
    read_qual = fq_handle.readline().strip()

    return read_header, read_seq, read_plus, read_qual

# Generate the index list and the index files dictionary
# Dictionary file format -> {index: [file_handle_r1, file_handle_r2]}
idx_list, idx_files_dict = create_index_files()

# Variable to keep track of read-pairs with matches, hops and unknowns
unk_count, match_count, hopped_count = 0, 0, 0

# Dictionary to keep track of every pair indexes
match_pair_count_dict = {}
hop_pair_count_dict = {}

with gzip.open(read1, 'rt') as r1, gzip.open(index1, 'rt') as i1, \
    gzip.open(read2, 'rt') as r2, gzip.open(index2, 'rt') as i2:

    while True:
        # Store records for Read 1
        r1_header, r1_seq, r1_plus, r1_qual = read_fq(r1)
        if r1_header == "": break
        # Store records for Index 1
        i1_header, i1_seq, i1_plus, i1_qual = read_fq(i1)

        # Store records for Read 2
        r2_header, r2_seq, r2_plus, r2_qual = read_fq(r2)

        # Store records for Index 2
        i2_header, i2_seq, i2_plus, i2_qual = read_fq(i2)

        # Find the reverse complement of idx2
        i2_seq_rc = reverse_complement(i2_seq)

        # Modify the header 
        i1_i2 = i1_seq + '-' + i2_seq_rc
        new_r1_header = r1_header + " " + i1_i2
        new_r2_header = r2_header + " " + i1_i2

        # Search for unknowns first
        if search_for_unks(idx1=i1_seq, idx2=i2_seq_rc, idx1_qual=i1_qual, idx2_qual=i2_qual,\
                           idx_list=idx_list, th=qual_index):
            # Write out to unknown_r1
            idx_files_dict['unk'][0].write(f'{new_r1_header}\n{r1_seq}\n{r1_plus}\n{r1_qual}\n')

            # Write out to unknown_r2
            idx_files_dict['unk'][1].write(f'{new_r2_header}\n{r2_seq}\n{r2_plus}\n{r2_qual}\n')

            unk_count += 1
            continue

        if search_for_hops(i1_seq, i2_seq_rc):
            # Write out to hop_r1
            idx_files_dict['hop'][0].write(f'{new_r1_header}\n{r1_seq}\n{r1_plus}\n{r1_qual}\n')

            # Write out to hop_r2
            idx_files_dict['hop'][1].write(f'{new_r2_header}\n{r2_seq}\n{r2_plus}\n{r2_qual}\n')
            hopped_count += 1
            # Add read-pair to dict
            if hop_pair_count_dict.get(i1_i2) == None:
                hop_pair_count_dict[i1_i2] = 1
            else:
                hop_pair_count_dict[i1_i2] += 1
            continue

        # Matching here:
        # Write out to match_r1
        idx_files_dict[i1_seq][0].write(f'{new_r1_header}\n{r1_seq}\n{r1_plus}\n{r1_qual}\n')

        # Write out to unknown_r2
        idx_files_dict[i1_seq][1].write(f'{new_r2_header}\n{r2_seq}\n{r2_plus}\n{r2_qual}\n')

        # Add read-pair to dict
        if match_pair_count_dict.get(i1_i2) == None:
            match_pair_count_dict[i1_i2] = 1
        else:
            match_pair_count_dict[i1_i2] += 1
        match_count += 1
        
# Write the report
total_reads = match_count + hopped_count + unk_count
report = open('report.txt', 'w')
report.write('REPORT SUMMARY (Percentage is Based over Total Reads)\n---------------------------------\n')
report.write(f'Matching Indexes:\t{match_count} records\t{match_count*100/total_reads:.2f}%\n\
Hopped Index:\t{hopped_count} records\t{hopped_count*100/total_reads}%\n\
Unkowns:\t{unk_count} records\t{unk_count*100/total_reads:.2f}%\n\
Total Reads: {total_reads}\n\
Qscore Cutoff: {qual_index}\n\n\
Matching Read Pairs Count (Percentage is Based over Properly Matched Pairs)\n---------------------------------\n')

for idx_pair in match_pair_count_dict:
    report.write(f'{idx_pair}:\t{match_pair_count_dict[idx_pair]}\t{match_pair_count_dict[idx_pair]*100/match_count:.2f}%\n')

report.write('\nHopped Read Pairs Count (Percentage is Based over Swapped Pairs)\n---------------------------------\n')

for idx_pair in hop_pair_count_dict:
    report.write(f'{idx_pair}:\t{hop_pair_count_dict[idx_pair]}\t{hop_pair_count_dict[idx_pair]*100/hopped_count:.2f}%\n')

report.close()
close_idx_files(idx_files_dict)