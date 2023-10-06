#!/usr/bin/env nextflow


process AVA{
    
    conda "bioconda::minimap2"
    cpus 8
    memory '4 GB'
    time 1.hour
        
    publishDir "${params.outdir}/${sample_name}", mode: 'copy'

    
    input: 
    val(sample_name)
    each path(splitted_reads)
    

    output:
    tuple val(sample_name), val(splitted_reads.baseName), path("${splitted_reads.baseName}.txt")

    
    script:
    """
    minimap2 -x ava-ont -D -t ${task.cpus} ${splitted_reads} ${splitted_reads} > ${splitted_reads.baseName}.txt
    """
}