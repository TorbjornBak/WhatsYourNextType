#!/usr/bin/env python3
import sys
import argparse
import subprocess
import os.path
from Bio import SeqIO
from Bio import pairwise2

class CustomParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)

def arguments():
    parser = CustomParser(description='Genotyper for HLA.')
    parser.add_argument('--fastq', type = str, required = True)
    parser.add_argument('--allelename', type = str, required= True)
    
   
    return parser.parse_args()

#Create a python script that iterates over a fastq file and splits reads according to some criterion based on their name

def pairwiseAlignment(fastqfile):
    for record in SeqIO.parse(fastqfile, "fastq"):
        


def main():
   args = arguments()
    
   

main()
