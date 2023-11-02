#!/usr/bin/env python3
import sys
import argparse
import subprocess
import os.path
import random
from Bio import SeqIO

def coverageCalculator(readdict, assemblylength):
    #print(readdict)
    reads = [len(sequence) for id, sequence in readdict.items()]
    #print(reads)
    return sum(reads) / assemblylength

def readDownSampler(readfile, assemblylength, coveragecutoff):
    #readdict = SeqIO.index(readfile, "fastq")
    readdict = {(read.id):(read.seq) for read in (SeqIO.parse(readfile, "fastq"))}

    coverage = coverageCalculator(readdict, assemblylength)
    print("Downsampling", readfile)
    print("Start coverage:", coverage)
    while coverage > coveragecutoff:
        for _ in range(int((coverage / coveragecutoff))):
            readdict.pop(random.choice(list(readdict.keys())))
        coverage = coverageCalculator(readdict, assemblylength)
    print("Final coverage:", coverage)
    return readdict

def writeDownsampledReads(readdict, outputfile,fastqfile):
    outfile = open(outputfile, 'w')
    
    for read in SeqIO.parse(fastqfile, "fastq"):
        if read.id in readdict:
            outfile.write(read.format('fastq'))


class CustomParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)


def arguments():
    parser = CustomParser(description='Genotyper for HLA.')
    parser.add_argument('--readfile', type = str)
    parser.add_argument('--outputfile', type = str)
    parser.add_argument('--assemblylength', type = int, default = 3000)
    parser.add_argument('--coveragecutoff', type = int, default = 100)
    
    return parser.parse_args()


def main():
    args = arguments()
    readdict = readDownSampler(args.readfile, args.assemblylength, args.coveragecutoff)
    writeDownsampledReads(readdict, args.outputfile, args.readfile)

main()