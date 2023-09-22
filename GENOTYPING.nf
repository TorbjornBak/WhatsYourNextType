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
    flye --nano-hq ${splitted_reads} --read-error 0.05 --threads ${task.cpus} --out-dir ${splitted_reads.baseName}
    """
}

process MINISAM{

    conda "bioconda::flye"
    
    cpus 8
    memory '4 GB'
    time 1.hour
        
    publishDir "${params.outdir}/${sample_name}", mode: 'copy'

    
    input: 
    tuple val(sample_name), val(allelename), path(assembly), path(splitted_reads)


    output:
    tuple val(sample_name), val(allelename), path(assembly), path("${allelename}_lr_mapping.bam")

    
    script:
    """
    minimap2 -ax map-ont -t ${task.cpus} ${assembly}/assembly.fasta ${splitted_reads} | samtools sort -@ 4 -m 1G > ${allelename}_lr_mapping.bam
    samtools index -@ 4 ${allelename}_lr_mapping.bam
    """


}

process HAPDUP{
    
    container = "mkolmogo/hapdup:0.2"
    cpus 8
    memory '4 GB'
    time 1.hour
        
    publishDir "${params.outdir}/${sample_name}", mode: 'copy'

    
    input: 
    tuple val(sample_name), val(allelename), path(assembly), path(bamfile)

    output:
    tuple val(sample_name), val(allelename), path("${allelename}/hapdup")

    
    script:
    """
    hapdup --assembly ${assembly}/assembly.fasta --bam ${bamfile} --out-dir ${allelename}/hapdup -t ${task.cpus} --rtype hifi
    """
}
