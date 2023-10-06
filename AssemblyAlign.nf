#!/usr/bin/env nextflow

process FLYEASSEMBLY{
    
    conda "bioconda::flye"
    cpus 8
    memory '4 GB'
    time 1.hour
        
    publishDir "${params.outdir}/${sample_name}", mode: 'copy'

    
    input: 
    val(sample_name)
    each path(splitted_reads)
    

    output:
    tuple val(sample_name), val(splitted_reads.baseName), path("${splitted_reads.baseName}"), path(splitted_reads)

    
    script:
    """
    flye --nano-hq ${splitted_reads} --read-error 0.05 --threads ${task.cpus} --out-dir ${splitted_reads.baseName} --min-overlap 2000 --genome-size 3500 --asm-coverage 100
    """
}

process MAFFT{
    
    conda "bioconda::mafft"
    cpus 8
    memory '4 GB'
    time 1.hour
        
    publishDir "${params.outdir}/${sample_name}", mode: 'copy'

    
    input: 
    tuple val(sample_name), val(allelename), path(assembly), path(splitted_reads)

    

    output:

    
    script:
    """
    mafft --6merpair --addfragments ${splitted_reads} ${assemly}/assembly.fasta > ${allelename}.mafft    
    """
}