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