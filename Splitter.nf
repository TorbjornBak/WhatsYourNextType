#!/usr/bin/env nextflow



process Splitter{

    cpus 8
    memory '4 GB'
    
    publishDir "${params.outdir}/${samplename}", mode: 'copy'

    
    input: 
    path(primerlist)
    path(fastqFile)
    val samplename

    output:
    path("Bins/*")
    val samplename
    
    script:
    """
    rm -rf -p Bins 
    mkdir Bins
    python PrimerSplitter.py $primerlist $fastqFile
    """
}