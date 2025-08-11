# 24/07/2025
## Data Exploration
4 fastq files are located in the Talapas directory: ``
`/projects/bgmp/shared/2017_sequencing/`

The directory is saved as the variable `dir` for easier calling.
```
dir=/projects/bgmp/shared/2017_sequencing/
```

Saved the following files in an array to later write a bash `for loop` for data exploration.
```
reads="
1294_S1_L008_R1_001.fastq.gz
1294_S1_L008_R2_001.fastq.gz
1294_S1_L008_R3_001.fastq.gz
1294_S1_L008_R4_001.fastq.gz
"
```

Bash pipeline used for data exploration:
```
for read in $reads; do echo ${read}; zcat ${dir}${read} | head -n 4; echo ""; done
```
The output yields:
```
1294_S1_L008_R1_001.fastq.gz
@K00337:83:HJKJNBBXX:8:1101:1265:1191 1:N:0:1
GNCTGGCATTCCCAGAGACATCAGTACCCAGTTGGTTCAGACAGTTCCTCTATTGGTTGACAAGGTCTTCATTTCTAGTGATATCAACACGGTGTCTACAA
+
A#A-<FJJJ<JJJJJJJJJJJJJJJJJFJJJJFFJJFJJJAJJJJ-AJJJJJJJFFJJJJJJFFA-7<AJJJFFAJJJJJF<F--JJJJJJF-A-F7JJJJ

1294_S1_L008_R2_001.fastq.gz
@K00337:83:HJKJNBBXX:8:1101:1265:1191 2:N:0:1
NCTTCGAC
+
#AA<FJJJ

1294_S1_L008_R3_001.fastq.gz
@K00337:83:HJKJNBBXX:8:1101:1265:1191 3:N:0:1
NTCGAAGA
+
#AAAAJJF

1294_S1_L008_R4_001.fastq.gz
@K00337:83:HJKJNBBXX:8:1101:1265:1191 4:N:0:1
NTTTTGATTTACCTTTCAGCCAATGAGAAGGCCGTTCATGCAGACTTTTTTAATGATTTTGAAGACCTTTTTGATGATGATGATGTCCAGTGAGGCCTCCC
+
#AAFAFJJ-----F---7-<FA-F<AFFA-JJJ77<FJFJFJJJJJJJJJJAFJFFAJJJJJJJJFJF7-AFFJJ7F7JFJJFJ7FFF--A<A7<-A-7--
```
Next, added commands to parse just the sequence and count the number of bases. Length of the read's sequence should be 108 then 8 for the barcode.
```bash
for read in $reads; do echo ${read}; zcat ${dir}${read} | head -n 4 | head -n 2 | tail -1 | wc -c; echo ""; done

1294_S1_L008_R1_001.fastq.gz
102

1294_S1_L008_R2_001.fastq.gz
9

1294_S1_L008_R3_001.fastq.gz
9

1294_S1_L008_R4_001.fastq.gz
102
```
Initial data exploration informs us that:
1294_S1_L008_R1_001.fastq.gz -> Read 1
1294_S1_L008_R2_001.fastq.gz -> Index 1
1294_S1_L008_R3_001.fastq.gz -> Index 2
1294_S1_L008_R4_001.fastq.gz  -> Read 2

### Determining pHRED Score
Earlier data exploration looking at the first record of each line contains the `>` `<` symbol which is not present in phred 64. Therefore, this is phred 33 encoding.

# Part 2
## Defining the Problem
The goal is to parse through 4 different fastq files to determine whether indexes of paired reads are matched, hopped or unknown. Another consideration taken into account is the index quality score as well. Our known index fq files are read 2 (the forward barcode) and read 3 (the reverse barcode). Read 3's sequence is the reverse complement of Read 2 as a result of the direction it was read prior to the bridge amplification.

The output of the script will give us 6 different fastq files. Firstly, a pair of fq files for reads that has properly matching index; secondly, a pair of fq files for reads with index-hopping observed; lastly, a pair of fq files for reads with either unknown indices or low quality index score.
### Pesudocode 
```python
def get_args():
	""" Parse the cmd-line args for 4 input files: R1, R2, R3, R4""

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

Bash Pipeline for creating unit tests:
```
i=1
for read in $reads; zcat ${dir}${read} | head -n 40 > test_r${i}; ((i++)); done
```
# 07/28/2025
Continued with data exploration for part 1.2. Since this is my first time using gzip in python, I decided to experiment with it to get a feel of what the syntax and output looks like.

Nothing out of the ordinary from the conventional approach of opening a file aside from having to include `gzip` with open.
```
with gzip.open(File_name, 'r') as fin:
	...
```
Tested with `1294_S1_L008_R1_001.fastq.gz ` using the following code:
```
read_file = '/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz'
with gzip.open(read_file, 'rb') as fin:
        print(fin.readline())
```
Terminal prints the following:
```
b'@K00337:83:HJKJNBBXX:8:1101:1265:1191 1:N:0:1\n'
```
So of course, I applied `strip()` and I got this:
```
b'@K00337:83:HJKJNBBXX:8:1101:1265:1191 1:N:0:1'
```
The `b` in front stands for binary. I wonder if I can slice individual character and turns out it outputs the ascii value of the character. This means I won't have to use my convert phred func and instead, directly subtract 33 later on as I write out the script.

## Writing the function `process_reads`
Wrote a function called `process_reads` to process each fastq files, using a dictionary of running sum strategy to compute the mean quality score distribution. This will prevent my slurm script from getting OOM kill.
```python
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
                    bp_scores = {i : 0 for i in range(len(seq))} #<- running sums for the histogram

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

```
I also added a `try` and `exception` statements just in case somewhere during the script, it encounters an error processing a file for safe measures.
## Writing the function `plot_distribution`
Wrote a function to plot the distribution of the outputted dictionary from the function `process_reads`. Added an additional customization to have grids as well.
```python
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
```
## Writing the Concurrent Features for Multiprocessing
This approach would allow my script utilize more than one when running. It starts by defining an executor handle to store the list of functions that has to be ran alongside the chosen parameters. Then using the output dictionary, it'll then get carried over to the `plot_distribution`.
```python
reads = [read1, read2, read3, read4]

with concurrent.futures.ProcessPoolExecutor() as executor:
    results = list(executor.map(process_read, reads))
    print(results)
    
    for i, bp_dt in enumerate(results, start=1):
        print(bp_dt)
        plot_distribution(bp_dt, label=i)
```
## Unit_test script
I wrote a simple bash script that creates unit_tests from the 4 input files by only taking the 10 records to ensure my code ran correctly without crashing.
```bash
#!/usr/bin/env bash

# create unit tests

dir=/projects/bgmp/shared/2017_sequencing/
reads="
1294_S1_L008_R1_001.fastq.gz
1294_S1_L008_R2_001.fastq.gz
1294_S1_L008_R3_001.fastq.gz
1294_S1_L008_R4_001.fastq.gz
"

i=1
for read in $reads;
    do
        zcat ${dir}${read} | head -n 40 > unit_test_read${i}.fq;
        gzip unit_test_read${i}.fq;
        ((i++))
    done
```
## SLURM script
Wrote the slurm script to ran my `qual_dist.py` script.
```bash

#!/usr/bin/env bash

#SBATCH --account=bgmp
#SBATCH --partition=bgmp
#SBATCH --cpus-per-task=5
#SBATCH --mem=16G
#SBATCH --nodes=1
#SBATCH --job-name=plot_hist.sh
#SBATCH --output=plot_hist_output.log
#SBATCH --error=plot_hist_error.log

dir=/projects/bgmp/shared/2017_sequencing/
read1=${dir}1294_S1_L008_R1_001.fastq.gz
read2=${dir}1294_S1_L008_R2_001.fastq.gz
read3=${dir}1294_S1_L008_R3_001.fastq.gz
read4=${dir}1294_S1_L008_R4_001.fastq.gz

# read1=unit_test_read1.fq.gz
# read2=unit_test_read2.fq.gz
# read3=unit_test_read3.fq.gz
# read4=unit_test_read4.fq.gz

mamba activate base
/usr/bin/time -v ./qual_dist.py -r1 $read1 -r2 $read2 -r3 $read3 -r4 $read4
```
Results are stored in the directory on talapas: `/projects/bgmp/dnhem/bioinfo/Bi622/Demultiplex/Assignment-the-first` with the following files:
```shell
(base) [dnhem@login4 Assignment-the-first]$ ls
Answers.md            plot.sh       qual_score_distribution_read1.png  README.md                    unit_test_read2.fq.gz
bioinfo.py            pseudo.txt    qual_score_distribution_read2.png  Unit_test_10_records.png     unit_test_read3.fq.gz
plot_hist_error.log   __pycache__   qual_score_distribution_read3.png  unit_test1_10_records.fq.gz  unit_test_read4.fq.gz
plot_hist_output.log  qual_dist.py  qual_score_distribution_read4.png  unit_test_read1.fq.gz        unit_test.sh
```
## Slurm Output
```

```
# Part 3
The index list is stored in the directory: `/projects/bgmp/shared/2017_sequencing/indexes_list`.
Initial data exploration of the index files outputs 24 unique indices with 5 different column headers:
```shell
(base) [dnhem@login4 Assignment-the-third]$ head $idx
sample  group   treatment       index   index sequence
1       2A      control B1      GTAGCGTA
2       2B      control A5      CGATCGAT
3       2B      control C1      GATCAAGG
4       2C      mbnl    B9      AACAGCGA
6       2D      mbnl    C9      TAGCCATG
7       2E      fox     C3      CGGTAATC
8       2F      fox     B3      CTCTGGAT
10      2G      both    C4      TACCGGAT
11      2H      both    A11     CTAGCTCA
```
This will later be used in my `get_args()` argument later as I will parse through the index file to obtain the list of indexes.

IMPORTANT: Read 3 has to be read as the reverse complement when comparing it with read 2!
# 07/30/2025
## Defining the get_args function
I defined two argument groups to make things look more presentable when the `help function` is called. I group them under `Fastq Files Stdin` and `Quality Threshold Control`. For the first argument group, we will take in all 4 of our read files followed with its index file. For the second argument group, we will take just the index threshold.
```python
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
```
## Defining Helper Functions
### Index Parser Functions
The function `create_index_files` was defined to parse through the index file input. The goal is to extract all the 24 relevant index and return the a string list of that index followed with an dictionary with the indexes as keys and 2 file handles for paired matches as values. 2 more dictionary keys are then created to store 2 more file handles each accounting for hopped indexes and unknowns from our read files.

The list of index is later used in a later helper function to check for known indexes and later help close all our open file handles with ease in the end.
```python
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
```
### Index Conditions Function
Two functions are defined below to ease our boolean logic when searching for unknown or hopped indexes. 

The first function `search_for_unks` check if the pair of indexes are classified as unknowns if they fall under at lease one of the following category:
- Contains `N`
- Not a valid index (Doesn't belong to the list of 24 known indexes)
- If one base falls under a quality score -> bad quality score index -> very good chance of incorrect sample
If all conditions are bypassed, the function will return false meaning our index are both of good quality and valid.

The second function `search_for_hops` simply compare if both index are equal and return True or False otherwise.
```python
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
```
### Reverse_Complement Function
This function will serve to find the reverse complement of read 3 when comparing it to read 2.

```python
def reverse_complement(seq:str):
    """
    Return the reverse complement of the sequence.
    """
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    comp = ''

    for char in seq:
        comp += complement[char]
    return comp[::-1]
```
### Read_Fq function
I wrote this function to easily parse every 4 lines of a fq file, taking advantage of the `readline` pointer (thank you C), so each time this is called in a loop, im storing new records each time to a variable for later usage.
```python
def read_fq(fq_handle):
    """
    Given a fq file handle, read the record.
    """
    read_header = fq_handle.readline().strip()
    read_seq = fq_handle.readline().strip()
    read_plus = fq_handle.readline().strip()
    read_qual = fq_handle.readline().strip()

    return read_header, read_seq, read_plus, read_qual

```
### Closing Files Function
```python
def close_idx_files(idx_dict:dict):
    """
    Close all Index Files Open.
    """
    for idx in idx_dict:
        print(idx)
        idx_dict[idx][0].close()
        idx_dict[idx][1].close()
```
## Performing the Algorithm
Started out with defining our global variables and initializing our data structures:
```python
# Allocating the arguments as a global variable int he script
args = get_args()
read1 = args.read1
index1 = args.index1
read2 = args.read2
index2 = args.index2
index_file = args.index_file
#qual_read = args.read_threshold
qual_index = int(args.index_threshold)

# Generate the index list and the index files dictionary
# Dictionary file format -> {index: [file_handle_r1, file_handle_r2]}
idx_list, idx_files_dict = create_index_files()

# Variable to keep track of read-pairs with matches, hops and unknowns
unk_count, match_count, hopped_count = 0, 0, 0

# Dictionary to keep track of every pair indexes
match_pair_count_dict = {}
hop_pair_count_dict = {}
```
The code below details how each record is being handled. We're storing records of each read file to their respective variables. `gzip` with `rt` is used to read bit files in text format (takes longer). I called `reverse_complement` on `i2_seq` to make sure it can be comparable with `i1_seq` later on.

Since our desired output files will have a modified header in which the pair of index will be added to it, the last lines show just that.
```python
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
```
The final section below details how each condition is called using our functions defined earlier. If whatever first condition is met (**unknown first then hop then finally match**), it'll be written to that corresponding file. Our count dictionary initialized in the beginning will take up the index pair as its key keep will be increment by 1 every time it encounters that specific index pair OR create a new key with value of 1 if encountering for the first time.
```python
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

        # Write out to match_r2
        idx_files_dict[i1_seq][1].write(f'{new_r2_header}\n{r2_seq}\n{r2_plus}\n{r2_qual}\n')

        # Add read-pair to dict
        if match_pair_count_dict.get(i1_i2) == None:
            match_pair_count_dict[i1_i2] = 1
        else:
            match_pair_count_dict[i1_i2] += 1
        match_count += 1
```
## Slurm Script
```bash
#!/usr/bin/env bash

#SBATCH --account=bgmp
#SBATCH --partition=bgmp
#SBATCH --cpus-per-task=5
#SBATCH --mem=16G
#SBATCH --nodes=1
#SBATCH --job-name=demultiplex.sh
#SBATCH --output=demultiplex_output.log
#SBATCH --error=demultiplex_error.log

dir=/projects/bgmp/shared/2017_sequencing/
read1=${dir}1294_S1_L008_R1_001.fastq.gz
read2=${dir}1294_S1_L008_R2_001.fastq.gz
read3=${dir}1294_S1_L008_R3_001.fastq.gz
read4=${dir}1294_S1_L008_R4_001.fastq.gz
index=${dir}indexes.txt

# read1=unit_test_read1.fq.gz
# read2=unit_test_read2.fq.gz
# read3=unit_test_read3.fq.gz
# read4=unit_test_read4.fq.gz


mamba activate base
/usr/bin/time -v ./demultiplex.py -r1 $read1 -r2 $read2 -r3 $read3 -r4 $read4 -idx $index -qi 30
```
## Writing the Report
Lastly, let's write up the report on the statistical distribution of our output keeping track of the matches, hops and unknowns. I'll be also taking into account the composition % of each as well.
```python
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
```
## Slurm Output for Cutoff 0
```
	Command being timed: "./demultiplex.py -r1 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz -r2 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz -r3 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz -r4 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz -idx /projects/bgmp/shared/2017_sequencing/indexes.txt -qi 0"
	User time (seconds): 3016.61
	System time (seconds): 47.11
	Percent of CPU this job got: 94%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 53:51.09
	Average shared text size (kbytes): 0
	Average unshared data size (kbytes): 0
	Average stack size (kbytes): 0
	Average total size (kbytes): 0
	Maximum resident set size (kbytes): 250908
	Average resident set size (kbytes): 0
	Major (requiring I/O) page faults: 0
	Minor (reclaiming a frame) page faults: 38256
	Voluntary context switches: 87550
	Involuntary context switches: 2357
	Swaps: 0
	File system inputs: 0
	File system outputs: 0
	Socket messages sent: 0
	Socket messages received: 0
	Signals delivered: 0
	Page size (bytes): 4096
	Exit status: 0

```
## Slurm Output for Cutoff 30
```
	Command being timed: "./demultiplex.py -r1 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz -r2 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz -r3 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz -r4 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz -idx /projects/bgmp/shared/2017_sequencing/indexes.txt -qi 30"
	User time (seconds): 2915.20
	System time (seconds): 34.89
	Percent of CPU this job got: 95%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 51:36.73
	Average shared text size (kbytes): 0
	Average unshared data size (kbytes): 0
	Average stack size (kbytes): 0
	Average total size (kbytes): 0
	Maximum resident set size (kbytes): 251768
	Average resident set size (kbytes): 0
	Major (requiring I/O) page faults: 0
	Minor (reclaiming a frame) page faults: 38333
	Voluntary context switches: 48331
	Involuntary context switches: 918
	Swaps: 0
	File system inputs: 0
	File system outputs: 0
	Socket messages sent: 0
	Socket messages received: 0
	Signals delivered: 0
	Page size (bytes): 4096
	Exit status: 0

```

I ran the code twice one for no cutoff and one for cutoff=30.
## Report with no cutoff:
```
REPORT SUMMARY (Percentage is Based over Total Reads)
---------------------------------
Matching Indexes:	331755033 records	91.33%
Hopped Index:	707740 records	0.19483726398807136%
Unkowns:	30783962 records	8.47%
Total Reads: 363246735
Qscore Cutoff: 0
```
## Report with Cutoff
```
REPORT SUMMARY (Percentage is Based over Total Reads)
---------------------------------
Matching Indexes:	226715602 records	62.41%
Hopped Index:	330975 records	0.09111575359376596%
Unkowns:	136200158 records	37.50%
Total Reads: 363246735
Qscore Cutoff: 30
```