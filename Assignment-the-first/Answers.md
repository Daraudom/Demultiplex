# Assignment the First

## Part 1
1. Be sure to upload your Python script. Provide a link to it here:

| File name | label | Read length | Phred encoding |
|---|---|---|---|
| 1294_S1_L008_R1_001.fastq.gz |  |  |  |
| 1294_S1_L008_R2_001.fastq.gz |  |  |  |
| 1294_S1_L008_R3_001.fastq.gz |  |  |  |
| 1294_S1_L008_R4_001.fastq.gz |  |  |  |

2. Per-base NT distribution
    1. Use markdown to insert your 4 histograms here.
    2. **YOUR ANSWER HERE**
    3. **YOUR ANSWER HERE**
    
## Part 2
1. Define the problem

The goal is to parse through 4 different fastq files to determine whether indexes of paired reads are matched, hopped or unknown. Another consideration taken into account is the index quality score as well. Our known index fq files are read 2 (the forward barcode) and read 3 (the reverse barcode). Read 3's sequence is the reverse complement of Read 2 as a result of the direction it was read prior to the bridge amplification.

2. Describe output

The output of the script will give us 6 different fastq files. Firstly, a pair of fq files for reads that has properly matching index; secondly, a pair of fq files for reads with index-hopping observed; lastly, a pair of fq files for reads with either unknown indices or low quality index score.

Additionally, there will be a generated report txt file that outlines the counts for matching index, index-hopping and unknowns.

3. Upload your [4 input FASTQ files](../TEST-input_FASTQ) and your [>=6 expected output FASTQ files](../TEST-output_FASTQ).
4. Pseudocode
```python
r1, r2, r3, r4 <- open 4 files to read
match_r1, match_r2, hop_r1, hop_r2, unk_r1, unk_r2 <- open 4 files to write

List_of_idx = [...]
counter = Dictionary with key:Value pairs of {match:0, hop:0, unk:0}
idx_pair_dict = Dictionary with key:Value pairs of {idx:0} <- tracks count of unique index pairs

while true:
	Read records of  r1, r2, r3, r4 <- track header, seq, +, qual as vars for each
	r1_head, r1_seq, r1_pls, r1_qual <- r1_records
	r2_head, r2_seq, r2_pls, r2_qual <- r2_records
	r3_head, r3_seq, r3_pls, r3_qual <- r3_records
	r4_head, r4_seq, r4_pls, r4_qual <- r4_records
	
	Break while loop once EOF is hit
	
	Find the reverse complement of r3
	r3_rc = reverse_complement(r3_seq)
	
	Convert phred score of r2 and r3 into integer:
	r2_idx_scr = convert_phred(r2); r3_idx_scr = convert_phred(r3_rc)
	# Both idx_scr will be a list of ints
	
	Update headers of index-pairs:
	idx_pair = r2_seq + '-' + r3_seq
	updated_r1_head = r1_head + idx_pair
	updated_r4_head = r4_head + idx_pair
	
	Check if idx_pair exist in idx_pair_dict
		-> create a key with value 1 if not existent
		-> increment existing key if it exists
	
	Assess the quality of r2 and r3 index score
	if r2 or r3 contains N OR missing from list of index OR < qual_threshold:
		unk_r1.write(r1_records); unk_r2.write(r4_records)
		-> incrememnt counter[unk]
	elif r2_seq not equal to r3_rc:
		hop_r1.write(r1_records); hop_r2.write(r4_records)
		-> incrememnt counter[hop]
	else:
		match_r1.write(r1_records); match_r2.write(r2_records)
		-> incrememnt counter[match]
	
	
```
5. High level functions. For each function, be sure to include:
    1. Description/doc string
    2. Function headers (name and parameters)
    3. Test examples for individual functions
    4. Return statement
```python
def get_args():
	'''Parse the cmd-line args for 4 input files: R1, R2, R3, R4'''

Input: 4 fq files: r1, r2, r3, 4
Expected Output: Parser.args()

def convert_phred(letter: str) -> int:
    '''Takes a single ASCII character (string) encoded in Phred+33 and
    returns the quality score value as an integer.'''
    return qscore
    
Input: I
Expected output: 40

def reverse_complement(seq: str) -> str:
	'''Takes a DNA sequence and returns the reverse complement.'''

Input: ACTG
Expected output: TGAC
```