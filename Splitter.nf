#!/usr/bin/env nextflow



process SPLITTER{

    conda = "bioconda::magicblast conda-forge::biopython"
    cpus 4
    memory '4 GB'
    publishDir "${params.outdir}/${samplename}", mode: 'copy'

    
    input: 
    //path(primerlist)
    path(fastqFile)
    val samplename

    output:
    tuple val(samplename), path("Bins/*_bin.fastq")
    
    
    script:
    """
    rm -rf -f Bins 
    mkdir Bins
    python3 ${projectDir}/PrimerSplitter.py ${projectDir}/${params.primerlist} ${fastqFile} \$PWD
    magicblast -query Bins/excess2_bin.fastq -db ${projectDir}/${params.blastdb} -out Bins/excess_bin_blast.txt -num_threads ${task.cpus} -gapopen 20 -gapextend 20 -infmt fastq -outfmt tabular
    python3 ${projectDir}/BlastSplitter.py --blastfile Bins/excess_bin_blast.txt --allelelist ${projectDir}/Allelelist.txt --hlagnom ${projectDir}/hla_nom_g.txt --fastq Bins/excess2_bin.fastq --oA Bins/HLAA_bin.fastq --oB Bins/HLAB_bin.fastq --oC Bins/HLAC_bin.fastq --oDRB1 Bins/DRB1_bin.fastq --oDQA1 Bins/DQA1_bin.fastq --oDQB1 Bins/DQB1_bin.fastq --oDPB1 Bins/DPB1_bin.fastq
    rm Bins/excess_bin.fastq -f
    rm Bins/excess2_bin.fastq -f
    """
}

