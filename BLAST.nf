#!/usr/bin/env nextflow

process BLASTN{
    
    conda "bioconda::blast"
    cpus 8
    memory '4 GB'
    time 1.hour
        
    publishDir "${params.outdir}/${sample_name}", mode: 'copy'

    
    input: 
    tuple val(sample_name), val(allelename), path(assembly), path(splitted_reads)
    

    output:
    tuple val(sample_name), val(allelename), path("${allelelname}_blastresults.txt")
    
    script:
    """
    blastn -query ${assembly} -db ${params.blastdb} -out ${allelelname}_blastresults.txt -num_alignments 3 -gapopen 3 -gapextend 2 -penalty 3 -reward 1 -word_size 11 -num_threads ${task.cpus}
    """
}