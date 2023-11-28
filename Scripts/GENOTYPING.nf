#!/usr/bin/env nextflow

process ASSEMBLY{
    errorStrategy {task.attempt < 4 ? 'retry' : 'ignore'}
    conda "bioconda::flye"
    time 1.hour
    maxRetries 4
    tag "${sample_name}:${splitted_reads.baseName}"
        
    // publishDir "${params.outdir}/${sample_name}/flye_assembly", mode: 'copy'

    
    input: 
    tuple val(sample_name), path(splitted_reads), path(expectedgenomesize)
    


    output:
    tuple val(sample_name), val(splitted_reads.baseName), path("${splitted_reads.baseName}"), path(splitted_reads), path("${splitted_reads.baseName}/assembly.fasta")
    

    
    script:
    if (task.attempt == 1) {
    """
    flye --nano-hq ${splitted_reads} --read-error ${params.readerror} --min-overlap 2000 --threads ${task.cpus} --out-dir ${splitted_reads.baseName} --genome-size ${expectedgenomesize.baseName} --iterations 3 --asm-coverage ${params.coverage}
    """
    }
    else if ((task.attempt == 2)) {
    """
    flye --nano-hq ${splitted_reads} --read-error ${params.readerror} --min-overlap 1800 --threads ${task.cpus} --out-dir ${splitted_reads.baseName} --genome-size ${expectedgenomesize.baseName} --iterations 3 --asm-coverage ${params.coverage}
    """
    }
    else if ((task.attempt == 3)) {
    """
    flye --nano-hq ${splitted_reads} --read-error ${params.readerror} --min-overlap 1500 --threads ${task.cpus} --out-dir ${splitted_reads.baseName} --genome-size ${expectedgenomesize.baseName} --iterations 3 --asm-coverage ${params.coverage}
    """
    }
    else{
    """
    flye --nano-hq ${splitted_reads} --read-error ${params.readerror} --min-overlap 1001 --threads ${task.cpus} --out-dir ${splitted_reads.baseName} --genome-size ${expectedgenomesize.baseName} --iterations 3 --asm-coverage ${params.coverage}
    """
    }


}


process MINISAM{

    conda "bioconda::minimap2 bioconda::samtools"
    
    cpus 8
    memory '4 GB'
    time 1.hour
    tag "${sample_name}:${allelename}"

    //publishDir "${params.outdir}/${sample_name}/minisam", mode: 'copy'

    
    input: 
    tuple val(sample_name), val(allelename), path(assembly), path(splitted_reads), val(contigs)
   

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


process FILTERNOTMAPPEDREADS{

    conda "bioconda::minimap2 bioconda::samtools"
    
    cpus 8
    memory '4 GB'
    time 1.hour
        
    publishDir "${params.outdir}/${sample_name}/minisam", mode: 'copy'

    
    input: 
    tuple val(sample_name), val(allelename), path(assembly), path(splitted_reads), val(contigs)
   

    output:
    tuple val(sample_name), val(allelename), path(assembly), path("${allelename}_lr_mapping.bam"), path("${allelename}_lr_mapping.bam.bai"), path(splitted_reads)

    
    script:
    """
    minimap2 -ax map-ont ${assembly}/assembly.fasta ${splitted_reads} | samtools fastq -n -f 4 - > unassembled_reads.fastq.gz
    """

}


process PHASING{
    errorStrategy 'retry'
    container = "mkolmogo/hapdup:0.12"
    cpus 4
    memory '12 GB'
    time 1.hour
    maxRetries 3
    tag "${sample_name}:${allelename}"

    
    publishDir "${params.outdir}/${sample_name}/hapdup/", mode: 'copy'

    
    input: 
    tuple val(sample_name), val(allelename), path(assembly), path(bamfile), path(indexfile), path(splitted_reads)

    output:
    tuple val(sample_name), val(allelename), path("${allelename}_hapdup/${allelename}_hapdup_dual_*.fasta"), path("${allelename}_hapdup/")
    tuple val(sample_name), val(allelename), path("${allelename}_hapdup/")


    script:
    if (task.attempt == 1) {
    """
    hapdup --assembly ${assembly}/assembly.fasta --bam ${bamfile} --bam-index ${indexfile} --out-dir ${allelename}_hapdup --rtype hifi -t ${task.cpus} --min-aligned-length 1200 --max-read-error 0.10 --overwrite

    cp ${allelename}_hapdup/hapdup_dual_1.fasta ${allelename}_hapdup/${allelename}_hapdup_dual_1.fasta
    cp ${allelename}_hapdup/hapdup_dual_2.fasta ${allelename}_hapdup/${allelename}_hapdup_dual_2.fasta
    cp ${allelename}_hapdup/margin/margin.log ${allelename}_hapdup/margin/${allelename}_margin.log
    """
    }
    else {
    """
    mkdir ${allelename}_hapdup -p
    cp ${assembly}/assembly.fasta ${allelename}_hapdup/${allelename}_hapdup_dual_nophase.fasta
    """
    }
}

process EXTRACTMARGIN {

    // publishDir "${params.outdir}/${sample_name}/", mode: 'copy'

    input:
    tuple val(sample_name), val(allelename), path(hapdupdir)

    output:
    tuple val(sample_name), path("${sample_name}_haplotypedistribution.csv")
    
    script:
    """
    for var in ${hapdupdir}
    do
    haplotypes=\$(grep Separated \$var/margin/*_margin.log | cut -d " " -f8,10,13)
    echo \$var, \$haplotypes >> ${sample_name}_haplotypedistribution.csv
    done
    """
}


process UNPHASED {
    tag "${sample_name}:${allelename}"

    input:
    tuple val(sample_name), val(allelename), path(assembly), path(fastq), val(contigs)

    output:
    tuple val(sample_name), val(allelename), path("${allelename}_assembly.fasta"), path(assembly)
    
    script:
    """

    
    cp ${assembly}/assembly.fasta ${allelename}_assembly.fasta
    """
}

