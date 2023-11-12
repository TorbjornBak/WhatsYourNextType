#!/usr/bin/env python3
import sys
import argparse
import subprocess
import os.path
import random
from Bio import SeqIO
import multiprocessing


def coverageCalculator(readdict, assemblylength):
    #print(readdict)
    reads = [len(sequence) for id, sequence in readdict.items()]
    #print(reads)
    return sum(reads) / assemblylength

def process_chunk(chunk, length_cutoff):
    filtered_reads = {(read.id): (read.seq) for read in chunk if len(read.seq) > length_cutoff}
    return filtered_reads

def readDownSampler(readfile, assemblylength, coveragecutoff, lengthcutoff, threads = multiprocessing.cpu_count()):


    # Define the number of processes to use
    num_processes = threads

    # Read the FASTQ file in chunks
    records = list(SeqIO.parse(readfile, "fastq"))
    chunk_size = len(records) // num_processes
    chunks = [records[i:i + chunk_size] for i in range(0, len(records), chunk_size)]

    # Create a multiprocessing Pool
    with multiprocessing.Pool(processes=num_processes) as pool:
        # Use multiprocessing to filter reads in parallel
        results = pool.starmap(process_chunk, [(chunk, lengthcutoff) for chunk in chunks])

    # Combine results from different processes
    readdict = {}
    for result in results:
        readdict.update(result)

    print("Filtered reads:", readdict)
    #Creating a dictionary containing ids and sequences of the reads
    #readdict = {(read.id):(read.seq) for read in (SeqIO.parse(readfile, "fastq")) if len(read.seq) > lengthcutoff}

    coverage = coverageCalculator(readdict, assemblylength)
    
    print("Downsampling", readfile)
    print("Assumed fragment length:", assemblylength)
    print("Start coverage:", coverage)
    while coverage > coveragecutoff:
        pops = random.sample(list(readdict.items()), k = int(coveragecutoff / coverage * len(list(readdict.items()))))

        readdict = {id:sequence for id, sequence in pops}
        coverage = coverageCalculator(readdict, assemblylength)
    print("Final coverage:", coverage)

    return readdict

def writeDownsampledReads(readdict, outputfile,fastqfile):
    outfile = open(outputfile, 'w')
    
    for read in SeqIO.parse(fastqfile, "fastq"):
        if read.id in readdict:
            outfile.write(read.format('fastq'))

def findFragmentLength(fragmentlengthfile, allelename):
    allele = allelename.split('_')[0]
    with open(fragmentlengthfile,'r') as file:
        for line in file:
            if line.startswith(allele):
                file = open(f'{line.split()[1]}.{allelename}','w')
                file.close()
                return int(line.split()[1])


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
    parser.add_argument('--fragmentlength', type = str, default = None)
    parser.add_argument('--allele', type = str, default = None)
    parser.add_argument('--readsizecutoff', type = int, default = 2000)
    parser.add_argument('--threads', type = int, default = multiprocessing.cpu_count())
    
    return parser.parse_args()


def main():
    args = arguments()
    if args.allele is not None:
        assemblylength = findFragmentLength(args.fragmentlength, args.allele)
    else:
        assemblylength = 3200

    readdict = readDownSampler(args.readfile, assemblylength, args.coveragecutoff, args.readsizecutoff, args.threads)
    writeDownsampledReads(readdict, args.outputfile, args.readfile)

main()