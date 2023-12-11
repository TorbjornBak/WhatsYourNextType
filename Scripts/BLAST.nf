#!/usr/bin/env nextflow

process BLASTN {
    
    //conda "bioconda::blast"
    cpus 8
    memory '4 GB'
    time 1.hour
    tag "${sample_name}:${allelename}"
        
    //publishDir "${params.outdir}/${sample_name}", mode: 'copy'

    
    input: 
    tuple val(sample_name), val(allelename), path(assembly), path(assemblyfolder)
    

    output:
    tuple val(sample_name), val(allelename), path("${assembly.baseName}_blastresults.txt")
    
    //script:
    //"""
    //blastn -query ${assembly}/assembly.fasta -db ${projectDir}/${params.blastdb} -out ${allelename}_blastresults.txt -num_alignments 3 -gapopen 3 -gapextend 2 -penalty -3 -reward 2 -word_size 28 -num_threads ${task.cpus}
    //"""
    script:
    """
    blastn -task megablast -query ${assembly} -db ${projectDir}/${params.blastdb} -out ${assembly.baseName}_blastresults.txt -num_alignments 3 -num_threads ${task.cpus} -gapopen 0 -gapextend 0
    """

}

process CATBLAST {

    cpus 1
    memory '4 GB'
    time 1.hour
    tag "${sample_name}:${allelename}"
    
    publishDir "${params.outdir}/${sample_name}/BlastResults", mode: 'copy'
      
    input: 
    tuple val(sample_name), val(allelename), path(blastresults)

    output:
    tuple val(sample_name), path("${sample_name}_Blastresults.txt")
    
    script:
    """
    cat *_blastresults.txt > ${sample_name}_Blastresults.txt
    """
}

process HLAGENOTYPER {
    debug true
    memory '4 GB'
    time 1.hour
    tag "${sample_name}"
    publishDir "${params.outdir}/${sample_name}", mode: 'copy'
    conda "pandas"

    
    input: 
    tuple val(sample_name), path(blastresults), path(haplotypedist)
    

    output:
    tuple val(sample_name), path ("${sample_name}_HLA_type.tsv")
    
    script:
    """

    python3 ${projectDir}/Scripts/HLA_genotyper.py --blastfile ${blastresults} --hlagen ${projectDir}/${params.hlaGfile} --output ${sample_name}_HLA_type.tsv --marginLog ${haplotypedist}
    """

}
