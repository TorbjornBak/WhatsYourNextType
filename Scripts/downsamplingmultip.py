#!/usr/bin/env python3
import sys
import argparse
import subprocess
import os.path
import random
from Bio import SeqIO
import multiprocessing
import gzip


def coverageCalculator(readdict, assemblylength):
    #print(readdict)
    reads = [len(sequence) for id, sequence in readdict.items()]
    #print(reads)
    return sum(reads) / assemblylength

def process_chunk(chunk, lowercutoff,uppercutoff):
    filtered_reads = {(read.id): (read.seq) for read in chunk if len(read.seq) > lowercutoff and len(read.seq) < uppercutoff}
    return filtered_reads

def is_gz_file(filepath):
    with open(filepath, 'rb') as test_f:
        return test_f.read(2) == b'\x1f\x8b'


def readDownSampler(readfile, assemblylength, coveragecutoff, lowercutoff, uppercutoff, threads = multiprocessing.cpu_count()): 

    # Define the number of processes to use
    num_processes = threads

    # Read the FASTQ file in chunks
    if is_gz_file(readfile):
    
        with gzip.open(readfile, "rt") as file:
            records = list(SeqIO.parse(file, "fastq"))
        if len(records) == 0:
            sys.exit(1)

    else:
        records = list(SeqIO.parse(readfile, "fastq"))
    
    chunk_size = len(records) // num_processes
    chunks = [records[i:i + chunk_size] for i in range(0, len(records), chunk_size)]

    # Create a multiprocessing Pool
    with multiprocessing.Pool(processes=num_processes) as pool:
        # Use multiprocessing to filter reads in parallel
        results = pool.starmap(process_chunk, [(chunk, lowercutoff, uppercutoff) for chunk in chunks])

    # Combine results from different processes
    readdict = {}
    for result in results:
        readdict.update(result)

    #print("Filtered reads:", readdict)
    #Creating a dictionary containing ids and sequences of the reads
    #readdict = {(read.id):(read.seq) for read in (SeqIO.parse(readfile, "fastq")) if len(read.seq) > lengthcutoff}

    coverage = coverageCalculator(readdict, assemblylength)
    
    print(f"Downsampling {readfile} to {coveragecutoff}.")
    print(f"Assumed fragment length: {assemblylength}")
    
    print(f"Start coverage: {int(coverage)}")
    while coverage > coveragecutoff:
        pops = random.sample(list(readdict.items()), k = int(coveragecutoff / coverage * len(list(readdict.items()))))

        readdict = {id:sequence for id, sequence in pops}
        coverage = coverageCalculator(readdict, assemblylength)
    print(f"Final coverage: {int(coverage)}")

    return readdict

def writeDownsampledReads(readdict, outputfile, fastqfile):
    
    outfile = open(outputfile, 'w')

    if is_gz_file(fastqfile):
        with gzip.open(fastqfile, "rt") as file:
            for read in SeqIO.parse(file, "fastq"):
                if read.id in readdict:
                    outfile.write(read.format('fastq'))
    else:
        for read in SeqIO.parse(fastqfile, "fastq"):
            if read.id in readdict:
                outfile.write(read.format('fastq'))
        
    return
    

def findFragmentLength(fragmentlengthfile, allelename):
    allele = allelename.split('_')[0]
    with open(fragmentlengthfile,'r') as file:
        for line in file:
            if line.startswith(allele):
                file = open(f'{line.split()[1]}.{allelename}','w')
                file.close()
                return int(line.split()[1])

def getDynamicCoverage(allelename,coverage):
    allelename = allelename.split("_")[0] #Removing the bin in DRB1_bin
    coverageDict = {"HLAA":1,"HLAB":1,"HLAC":1,"DRB1":1,"DQA1":1,"DQB1":1,"DPB1":1}
    return int(coverageDict[allelename]*coverage)

class CustomParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)


def arguments():
    parser = CustomParser(description='Genotyper for HLA.')
    parser.add_argument('--readfile', type = str)
    parser.add_argument('--outputfile', type = str)
    parser.add_argument('--coveragecutoff', type = int, default = 100)
    parser.add_argument('--coveragedynamic', action = 'store_true')
    parser.add_argument('--fragmentlength', type = str, default = None)
    parser.add_argument('--allele', type = str, default = None)
    parser.add_argument('--lowercutoff', type = int, default = 2000)
    parser.add_argument('--uppercutoff', type = int, default = 3700)
    parser.add_argument('--threads', type = int, default = multiprocessing.cpu_count())
    
    return parser.parse_args()


def main():
    args = arguments()
    
    if args.allele is not None:
        assemblylength = findFragmentLength(args.fragmentlength, args.allele)
    else:
        assemblylength = 3200
    
    if args.coveragedynamic:
        coverage = getDynamicCoverage(args.allele,args.coveragecutoff)
    else:
        coverage = args.coveragecutoff

    readdict = readDownSampler(args.readfile, assemblylength, coverage, args.lowercutoff, args.uppercutoff, args.threads)
    writeDownsampledReads(readdict, args.outputfile, args.readfile)

main()
