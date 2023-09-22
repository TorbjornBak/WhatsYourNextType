#!/usr/bin/env nextflow

process FLYEASSEMBLY{
    
    conda "bioconda::flye"
    cpus 8
    memory '4 GB'
    time 1.hour
        
    publishDir "${params.outdir}/${sample_name}", mode: 'copy'

    
    input: 
    val(sample_name)
    path(splitted_reads)

    output:
    val(sample_name)
    path("Assembly/${splitted_reads.baseName()}")

    
    script:
    """
    flye --nano-hq ${splitted_reads} --read-error 0.05 --threads ${task.cpus} --out-dir Assembly/${splitted_reads.baseName()}
    """
}