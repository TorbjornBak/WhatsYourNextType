#!/bin/bash

# Check for the presence of required tools
command -v flye >/dev/null 2>&1 || { echo >&2 "Flye is not installed. Please install it before running this script."; exit 1; }

# Check if the input file is provided as an argument
if [ $# -ne 1 ]; then
    echo "Usage: $0 <input.fastq>"
    exit 1
fi

# Input FASTQ file
input_fastq="$1"

# Create a directory to store the assembly results
assembly_dir=assembly


# Run Flye assembly
flye --nano-hq "$input_fastq" --out-dir "$assembly_dir" --read-error 0.05 --threads 8 --min-overlap 2000 --genome-size 3500 --asm-coverage 100

blastn -task megablast -query "$assembly_dir"/assembly.fasta -db ../../HLAdatabase/HLA_db -out Test_blast_"$input_fastq".txt -num_alignments 3 -num_threads 4 -gapopen 20 -gapextend 20
