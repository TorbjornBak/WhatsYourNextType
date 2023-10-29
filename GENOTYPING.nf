#!/usr/bin/env nextflow

process FLYEASSEMBLY{
    //errorStrategy 'ignore'
    conda "bioconda::flye"
    cpus 8
    memory '4 GB'
    time 1.hour
        
    publishDir "${params.outdir}/${sample_name}", mode: 'copy'

    
    input: 
    tuple val(sample_name), path(splitted_reads)
    

    output:
    tuple val(sample_name), val(splitted_reads.baseName), path("${splitted_reads.baseName}"), path(splitted_reads)

    
    script:
    """
    flye --nano-hq ${splitted_reads} --read-error 0.05 --threads ${task.cpus} --out-dir ${splitted_reads.baseName} --min-overlap 2000 --genome-size 3000 --asm-coverage 800
    """
}



process SHASTA{
    
    conda "bioconda::shasta"
    cpus 8
    memory '4 GB'
    time 1.hour
        
    publishDir "${params.outdir}/${sample_name}", mode: 'copy'

    
    input: 
    tuple val(sample_name), path(splitted_reads)
    

    output:
    tuple val(sample_name), val(splitted_reads.baseName), path("${splitted_reads.baseName}"), path(splitted_reads)

    
    script:
    """
    shasta --config Nanopore-May2022 --input ${splitted_reads} --assemblyDirectory ${splitted_reads.baseName} --threads ${task.cpus} --Reads.minReadLength 2000 --Assembly.detangleMethod 2
    """
}



process FLYEASSEMBLY2{
    
    conda "bioconda::flye"
    cpus 8
    memory '4 GB'
    time 1.hour
        
    publishDir "${params.outdir}/${sample_name}", mode: 'copy'

    
    input: 
    tuple val(sample_name), val(allelename), path(splitted_reads), path(cluster1),path(cluster2)

    

    output:
    tuple val(sample_name), val(splitted_reads.baseName), path("${cluster1.baseName}"), path("${cluster2.baseName}"), path(splitted_reads)

    
    script:
    """
    flye --nano-hq ${splitted_reads} --read-error 0.05 --threads ${task.cpus} --out-dir ${cluster1.baseName} --min-overlap 2000 --genome-size 3500 --asm-coverage 100
    flye --nano-hq ${splitted_reads} --read-error 0.05 --threads ${task.cpus} --out-dir ${cluster2.baseName} --min-overlap 2000 --genome-size 3500 --asm-coverage 100

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
    tuple val(sample_name), val(allelename), path(assembly), path("${allelename}_lr_mapping.bam"), path("${allelename}_lr_mapping.bam.bai"), path(splitted_reads)

    
    script:
    """
    minimap2 -ax map-ont -t ${task.cpus} ${assembly}/assembly.fasta ${splitted_reads} | samtools sort -@ 4 -m 1G  > ${allelename}_lr_mapping_old.bam
    samtools addreplacerg -r '@RG\tID:${allelename}\tSM:${allelename}' -o ${allelename}_lr_mapping.bam ${allelename}_lr_mapping_old.bam
    samtools index -@ 4 ${allelename}_lr_mapping.bam
    samtools faidx ${assembly}/assembly.fasta 
    """
    // """
    // minimap2 -ax map-ont -t ${task.cpus} ${assembly}/assembly.fasta ${splitted_reads} | samtools sort -@ 4 -m 1G  > ${allelename}_lr_mapping.bam
    
    // samtools index -@ 4 ${allelename}_lr_mapping.bam
    
    // """


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
    errorStrategy 'retry'
    container = "mkolmogo/hapdup:0.12"
    cpus 4
    memory '12 GB'
    time 1.hour
    maxRetries 3

    
    publishDir "${params.outdir}/${sample_name}", mode: 'copy'

    
    input: 
    tuple val(sample_name), val(allelename), path(assembly), path(bamfile), path(indexfile), path(splitted_reads)

    output:
    tuple val(sample_name), val(allelename), path("${allelename}_hapdup/${allelename}_hapdup_dual_*.fasta"), path(splitted_reads)


    script:
    if (task.attempt == 1) {
    """
    hapdup --assembly ${assembly}/assembly.fasta --bam ${bamfile} --bam-index ${indexfile} --out-dir ${allelename}_hapdup --rtype hifi -t ${task.cpus} --min-aligned-length 1700 --max-read-error 0.07 --overwrite

    mv ${allelename}_hapdup/hapdup_dual_1.fasta ${allelename}_hapdup/${allelename}_hapdup_dual_1.fasta
    mv ${allelename}_hapdup/hapdup_dual_2.fasta ${allelename}_hapdup/${allelename}_hapdup_dual_2.fasta
    """
    }
    else {
    """
    mkdir ${allelename}_hapdup -p
    mv ${assembly}/assembly.fasta ${allelename}_hapdup/${allelename}_hapdup_dual_nophase.fasta
    """
    }
}

process WHATSHAP{
    //errorStrategy 'ignore'
    conda = "bioconda::whatshap"
    cpus 8
    memory '6 GB'
    time 1.hour
        
    publishDir "${params.outdir}/${sample_name}", mode: 'copy'

    
    input: 
    tuple val(sample_name), val(allelename), path(assembly), path(inputvcf), path(bamfile)

    output:
    

    script:
    """
    whatshap phase -o ${allelename}_phased.vcf --no-reference ${inputvcf} ${bamfile}
    """
}