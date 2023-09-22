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
    tuple val(sample_name), val(splitted_reads.baseName), path("Assembly/${splitted_reads.baseName}"), path(splitted_reads)

    
    script:
    """
    flye --nano-hq ${splitted_reads} --read-error 0.05 --threads ${task.cpus} --out-dir Assembly/${splitted_reads.baseName}
    """
}

process MINISAM{

    conda "bioconda::flye"
    
    cpus 8
    memory '4 GB'
    time 1.hour
        
    publishDir "${params.outdir}/${sample_name}", mode: 'copy'

    
    input: 
    tuple val(sample_name), val(allelename), path(assembly)


    output:
    tuple val(sample_name), val(allelename), path(assembly), path("${allelename}_lr_mapping.bam")

    
    script:
    """
    minimap2 -ax map-ont -t ${task.cpus} ${assembly}/assembly.fasta ${splitted_reads} | samtools sort -@ 4 -m 4G > ${allelename}_lr_mapping.bam
    samtools index -@ 4 ${allelename}_lr_mapping.bam
    """


}

process HAPDUP{
    
    container = "mkolmogo/hapdup"
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
