#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include {FLYEASSEMBLY,MINISAM,HAPDUP} from "./GENOTYPING.nf"
include {SPLITTER} from "./Splitter.nf"


// Main workflow script for the pipeline

workflow{

    Channel.fromPath(params.primerlist).set{primer_ch}
    Channel.fromPath(params.fastqfile).set{fastq_ch}
    
    SPLITREADS_ch = SPLITTER(primer_ch,fastq_ch,params.samplename)

    ASSEMBLY_ch = FLYEASSEMBLY(SPLITREADS_ch)
    MINISAM_ch = MINISAM(ASSEMBLY_ch)
    HAPDUP_ch = HAPDUP(MINISAM_ch)

    
}
