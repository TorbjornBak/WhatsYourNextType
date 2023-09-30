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

process MINISAM{

    conda "bioconda::minimap2 bioconda::samtools"
    
    cpus 8
    memory '4 GB'
    time 1.hour
        
    publishDir "${params.outdir}/${sample_name}", mode: 'copy'

    
    input: 
    tuple val(sample_name), val(allelename), path(assembly), path(splitted_reads)


    output:
    tuple val(sample_name), val(allelename), path(assembly), path("${allelename}_lr_mapping.bam"), path("${allelename}_lr_mapping.bam.bai")

    
    script:
    """
    minimap2 -ax map-ont -t ${task.cpus} ${assembly}/assembly.fasta ${splitted_reads} | samtools sort -@ 4 -m 1G  > ${allelename}_lr_mapping_old.bam
    samtools addreplacerg -r '@RG\tID:${allelename}\tSM:${allelename}' -o ${allelename}_lr_mapping.bam ${allelename}_lr_mapping_old.bam
    samtools index -@ 4 ${allelename}_lr_mapping.bam
    samtools faidx ${assembly}/assembly.fasta 
    """


}

process OCTOPUS{
    container "dancooke/octopus:latest"
    cpus 8
    memory '4 GB'
    time 1.hour
        
    publishDir "${params.outdir}/${sample_name}", mode: 'copy'

    
    input: 
    tuple val(sample_name), val(allelename), path(assembly), path(bamfile), path(indexfile)

    output:
    tuple val(sample_name), val(allelename), path(assembly), path("${allelename}.vcf.gz")

    
    script:
    """
    octopus -X 2G -B 6G --threads 8 --reads ${bamfile} --reference ${assembly}/assembly.fasta --disable-read-preprocessing -x 2 --dont-protect-reference-haplotype -o ${allelename}.vcf.gz
    """
}
process VCF{
    conda "vcftools"
    cpus 8
    memory '4 GB'
    time 1.hour
    publishDir "${params.outdir}/${sample_name}", mode: 'copy'

    input:
    tuple val(sample_name), val(allelename), path(assembly), path(VCFfile)

    output:
    tuple val(sample_name), val(allelename), path(assembly), path("${allelename}.consensus.fasta")

    script:
    """
    cat ${assembly}/assembly.fasta | vcf-consensus ${VCFfile} > ${allelename}.consensus.fasta
    """
}

process HAPDUP{
    container = "mkolmogo/hapdup:0.2"
    cpus 8
    memory '4 GB'
    time 1.hour
        
    publishDir "${params.outdir}/${sample_name}", mode: 'copy'

    
    input: 
    tuple val(sample_name), val(allelename), path(assembly), path(bamfile), path(indexfile)

    output:
    tuple val(sample_name), val(allelename), path("${allelename}/hapdup")

    
    script:
    """
    hapdup --assembly ${assembly}/assembly.fasta --bam ${bamfile} --out-dir ${allelename}/hapdup -t ${task.cpus}
    """
}
