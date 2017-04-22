#!/bin/bash
INPREFIX="/dsu0/ResearchData/174-12-26_READS/174-12-26-1-28280253/"
MATE1="${INPREFIX}174-12-26-1_S5_L001_R1_001.fastq.gz,${INPREFIX}174-12-26-1_S5_L002_R1_001.fastq.gz,${INPREFIX}174-12-26-1_S5_L003_R1_001.fastq.gz,${INPREFIX}174-12-26-1_S5_L004_R1_001.fastq.gz"
MATE2="${INPREFIX}174-12-26-1_S5_L001_R2_001.fastq.gz,${INPREFIX}174-12-26-1_S5_L002_R2_001.fastq.gz,${INPREFIX}174-12-26-1_S5_L003_R2_001.fastq.gz,${INPREFIX}174-12-26-1_S5_L004_R2_001.fastq.gz"
OUTFILE="174-12-26-1-28280253.bra.bowtie2.sam"
bowtie2 --fr --no-mixed --no-discordant --no-contain --no-overlap --no-unal -p 8 -x bra -1 $MATE1 -2 $MATE2 -S $OUTFILE

# OPTIONS EXPLAINED

# --fr	The upstream/downstream mate orientations for a valid paired-end 
# alignment against the forward reference strand. E.g., if --fr is specified and
# there is a candidate paired-end alignment where mate 1 appears upstream of the
# reverse complement of mate 2 and the fragment length constraints (-I and -X)
# are met, that alignment is valid.  Default: --fr (appropriate for Illumina's
# Paired-end Sequencing Assay).

# --no-mixed	By default, when bowtie2 cannot find a concordant or discordant
# alignment for a pair, it then tries to find alignments for the individual
# mates. This option disables that behavior.

# --no-discordant	By default, bowtie2 looks for discordant alignments if it
# cannot find any concordant alignments. A discordant alignment is an alignment
# where both mates align uniquely, but that does not satisfy the paired-end
# constraints (--fr/--rf/--ff, -I, -X). This option disables that behavior.

# --no-contain	If one mate alignment contains the other, consider that to be
# non-concordant. See also: Mates can overlap, contain or dovetail each other.
# Default: a mate can contain the other in a concordant alignment.

# --no-overlap	If one mate alignment overlaps the other at all, consider that
# to be non-concordant. See also: Mates can overlap, contain or dovetail each
# other. Default: mates can overlap in a concordant alignment.

# --no-unal	Suppress SAM records for reads that failed to align.

# -p/--threads NTHREADS	Launch NTHREADS parallel search threads (default: 1).
# Threads will run on separate processors/cores and synchronize when parsing
# reads and outputting alignments. Searching for alignments is highly parallel,
# and speedup is close to linear. Increasing -p increases Bowtie 2's memory
# footprint. E.g. when aligning to a human genome index, increasing -p from 1 to
# 8 increases the memory footprint by a few hundred megabytes. This option is
# only available if bowtie is linked with the pthreads library
# (i.e. if BOWTIE_PTHREADS=0 is not specified at build time).
