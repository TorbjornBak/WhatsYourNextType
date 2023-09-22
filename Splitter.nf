#!/usr/bin/env nextflow



process SPLITTER{

    cpus 8
    memory '4 GB'
    
    publishDir "${params.outdir}/${samplename}", mode: 'copy'

    
    input: 
    path(primerlist)
    path(fastqFile)
    val samplename

    output:
    val samplename
    path("Bins/*")
    
    
    script:
    """
    rm -rf -p Bins 
    mkdir Bins
    python ${projectDir}/PrimerSplitter.py $primerlist $fastqFile

    """
}