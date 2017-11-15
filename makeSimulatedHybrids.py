#!/usr/bin/env python

'''
* Daniel J. Shea
* F15M006C
* Niigata University, Dept. of Agriculture
* Laboratory of Plant Breeding

**2016-04-08:** This program will create simulated hybridized genomes for testing given the bol.fa and bra.fa fasta
reference sequences and a file with known syntenic regions for each of the chromosomes in the bol and bra genomes.

Copyright (c) 2016, Daniel J. Shea, Niigata University, Dept. of Agriculture, Laboratory of Plant Breeding
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the
following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following
   disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the
   following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote
   products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
'''

from Bio import SeqIO, SeqRecord
import random
import argparse
import sys

def parseSynteny(syntenyEntry):
    '''
    Split up the entry read in from the synteny file and convert coordinates from 1-based
    to 0-based python slice coordinates.
    '''
    # Split syntenyEntry into the following fields
    #Num Chromosome Region Strand MappedChromosome MappedRegion
    recnum, bChr, bRange, bStrand, iChr, iRange = syntenyEntry.split('\t')
    bStart, bStop = bRange.split('..')
    bStart = int(bStart)-1 # convert from 1-based to 0-based
    bStop  = int(bStop)    # keep the same since we use slices to extract syntenic blocks
    iStart, iStop = iRange.split('..')
    iStart = int(iStart)-1 # convert from 1-based to 0-based
    iStop  = int(iStop)    # keep the same since we use slices to extract syntenic blocks
    return (bChr, bStart, bStop, bStrand, iChr, iStart, iStop)


def makeIntrogressedGenome(outfilePrefix='IntrogressedGenome',
                           backgroundGenome='./simulate_genome/reference_genomes/bra.fa',
                           introgressGenome='./simulate_genome/reference_genomes/bol.fa',
                           syntenyFile='./simulate_genome/reference_genomes/bra_bol.synteny.txt',
                           seed=1234
                           ):
    '''
    Given two reference genome fasta files and the corresponding synteny file that lists
    all syntenic regions between the genomes, generate a genome that contains
    the replacement of one or more syntenic blocks with their corresponding regions.
    Write the new genome to an output file and record what changes were made to a separate
    file that can be used later to test if IntroMap got it right.
    '''
    with open(syntenyFile, 'r') as fh:
        fh.readline()  # Skips the header
        # Store synteny info into a list so we can randomly choose entries
        synteny = list()
        for line in fh:
            # Strip off the newline
            line = line.strip()
            synteny.append(line)

    outputLog = outfilePrefix + '.log'
    with open(outputLog, 'w') as fh:
        # Now, we select a random set of 5 to 10 syntenic blocks for the introgression
        # Although we want random data, we explicitly set a seed so the randomness is reproducible
        # We will record the seed used in the log file.
        random.seed(seed)
        fh.write('Seed set to {} for this run\n'.format(seed))

        # First determine how many introgressed regions to introduce
        numRegions = random.choice(range(5, 11))
        fh.write('Number of regions chosen was {}\n'.format(numRegions))

        # Store the chosen regions in a set so we do not have duplicates
        syntenySet = set()
        for i in range(numRegions):
            sChoice = random.choice(synteny)
            fh.write('Regions chosen was:\n{}\n'.format(sChoice))
            syntenySet.add(sChoice)

        # Parse the selections with parseSynteny() and store the returned tuples in a list
        syntenyBlocks = list()
        for s in syntenySet:
            syntenyTuple = parseSynteny(s)
            syntenyBlocks.append(syntenyTuple)

        # Load the genomes into dictionaries where the key is the seqid
        bSeqRecords = dict()
        for seqrec in SeqIO.parse(backgroundGenome, 'fasta'):
            bSeqRecords[seqrec.id] = seqrec
        fh.write('Loaded {} as background genome\n'.format(backgroundGenome))

        iSeqRecords = dict()
        for seqrec in SeqIO.parse(introgressGenome, 'fasta'):
            iSeqRecords[seqrec.id] = seqrec
        fh.write('Loaded {} as introgression genome\n'.format(introgressGenome))

        outputGenome = outfilePrefix + '.fasta'
        fh.write('Creating {}\n'.format(outputGenome))

        # Modify backgroundGenome in place using chosen syntenic regions from the introgressGenome
        # Perform modifications in order
        for s in syntenyBlocks:
            # Break out the tuple values so we can extract the slice from introgressGenome
            bChr, bStart, bStop, bStrand, iChr, iStart, iStop = s
            iSeq = iSeqRecords[iChr][bStart:bStop].seq
            # reverse complement the sequence if it maps back to the '-' strand in backgroundGenome
            if bStrand == '-':
                iSeq = iSeq.reverse_complement()
            # Perform the introgression
            bSeqRecords[bChr] = bSeqRecords[bChr][:bStart] + iSeq + bSeqRecords[bChr][bStop - 1:]
            fh.write('Introgression from {}({}..{}) into {} on strand {} at ({}..{})\n'.format(iChr, iStart + 1,
                                                                                               iStop, bChr,
                                                                                               bStrand, bStart + 1,
                                                                                               bStop))
        fh.write('Introgressions completed, writing out to file {}\n'.format(outputGenome))
        result = SeqIO.write(bSeqRecords.values(), outputGenome, 'fasta')
        fh.write('{} records written\n'.format(result))

# Here, we define the main function which will make use of argparse to provide a command line version of the
# Jupyter Notebook implementation of makeSimulatedHybrids. I highly recommend looking at the notebook if you're
# reading the code, it is more thoroughly documented.
def main():
    parser = argparse.ArgumentParser(description='makeSimulatedHybrids',
                                     epilog='If you use IntroMap, please consider citing our paper.')
    parser.add_argument('-b', '--background', type=str, required=True, dest='background',
                        help='The reference FASTA file of the recurrent parent. i.e. - The genomic background.')
    parser.add_argument('-i', '--introgression', type=str, required=True, dest='introgression',
                        help='The reference FASTA file of the donor parent. i.e. - The introgressed genome.')
    parser.add_argument('-f', '--synteny', type=str, required=True, dest='synteny',
                        help='The file that contains syntenic regions between the recurrent and donor parents.')
    parser.add_argument('-o', '--outprefix', type=str, required=True, dest='outprefix',
                        help='The output prefix for the resulting FASTA file.')
    parser.add_argument('-s', '--seed', type=int, required=True, dest='seed',
                        help='Seed for the psuedo-random number generator.')
    args = parser.parse_args()
    outfilePrefix = args.outprefix + '-' + str(args.seed)
    makeIntrogressedGenome(outfilePrefix=outfilePrefix,
                           backgroundGenome=args.background,
                           introgressGenome=args.introgression,
                           syntenyFile=args.synteny,
                           seed=args.seed
                           )
    sys.exit(0)

# Allow IntroMap to be run as a command line utility.
if __name__ == '__main__':
    main()
