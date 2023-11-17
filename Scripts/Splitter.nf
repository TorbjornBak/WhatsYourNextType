#!/usr/bin/env nextflow
process FIRSTDOWNSAMPLING  {
    conda "bioconda::biopython"
    cpus 4
    memory '4 GB'
    time 1.hour
    tag "${fastqFile.simpleName}"
    
    // publishDir "${params.outdir}/${fastqFile.baseName}", mode: 'copy'

    
    input: 
    path(fastqFile)
    

    output:
    tuple val("${fastqFile.simpleName}"), path("${fastqFile.simpleName}_sub.fastq")

    
    script:
    """
    python3 ${projectDir}/Scripts/downsamplingmultip.py --readfile ${fastqFile} --outputfile ${fastqFile.simpleName}_sub.fastq --coveragecutoff 20000 --readsizecutoff 1800 --threads ${task.cpus}
    """
}


process SPLITTER{
    errorStrategy = "ignore"
    conda = "bioconda::magicblast conda-forge::biopython"
    cpus 4
    memory '4 GB'
    tag "${sample_name}"

    // publishDir "${params.outdir}/${sample_name}", mode: 'copy'

    
    input: 
    tuple val(sample_name), path(fastqFile)
    

    output:
    tuple val(sample_name), path("Bins/*_bin.fastq")
    
    
    script:
    """
    mkdir Bins
    python3 ${projectDir}/Scripts/PrimerSplitter.py ${projectDir}/${params.primerlist} ${fastqFile} \$PWD
    magicblast -query Bins/excess2_bin.fastq -db ${projectDir}/${params.blastdb} -out Bins/excess_bin_blast.txt -num_threads ${task.cpus} -gapopen 20 -gapextend 20 -infmt fastq -outfmt tabular
    python3 ${projectDir}/Scripts/BlastSplitter.py --blastfile Bins/excess_bin_blast.txt --allelelist ${projectDir}/Data/Allelelist.txt --hlagnom ${projectDir}/${params.hlaGfile} --fastq Bins/excess2_bin.fastq --oA Bins/HLAA_bin.fastq --oB Bins/HLAB_bin.fastq --oC Bins/HLAC_bin.fastq --oDRB1 Bins/DRB1_bin.fastq --oDQA1 Bins/DQA1_bin.fastq --oDQB1 Bins/DQB1_bin.fastq --oDPB1 Bins/DPB1_bin.fastq
    rm Bins/excess_bin.fastq -f
    rm Bins/excess2_bin.fastq -f
    """
}

    