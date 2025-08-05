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