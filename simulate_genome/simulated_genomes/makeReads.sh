#!/bin/bash

# Program: wgsim (short read simulator)
# Version: 0.3.1-r13
# Contact: Heng Li <lh3@sanger.ac.uk>
# 
# Usage:   wgsim [options] <in.ref.fa> <out.read1.fq> <out.read2.fq>
# 
# Options: -e FLOAT      base error rate [0.020]
#          -d INT        outer distance between the two ends [500]
#          -s INT        standard deviation [50]
#          -N INT        number of read pairs [1000000]
#          -1 INT        length of the first read [70]
#          -2 INT        length of the second read [70]
#          -r FLOAT      rate of mutations [0.0010]
#          -R FLOAT      fraction of indels [0.15]
#          -X FLOAT      probability an indel is extended [0.30]
#          -S INT        seed for random generator [-1]
#          -A FLOAT      disgard if the fraction of ambiguous bases higher than FLOAT [0.05]
#          -h            haplotype mode

wgsim=/dsu0/LocalSoftware/wgsim/wgsim
refs="20160409AC-1234.fasta 20160409AC-2468.fasta 20160409AC-36912.fasta"
for ref in $refs; do
	outfile=$(echo $ref | cut -d. -f1)
	outfile1="${outfile}_S5_L001_R1_001.fastq"
	outfile2="${outfile}_S5_L001_R2_001.fastq"
	$wgsim -1 75 -2 75 -N 50000000 -S 1234 $ref $outfile1 $outfile2
done
