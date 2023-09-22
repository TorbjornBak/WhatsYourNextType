#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include {FLYEASSEMBLY} from "./GENOTYPING.nf"
include {SPLITTER} from "./Splitter.nf"


// Main workflow script for the pipeline

workflow{

    SPLITREADS_ch = SPLITTER(params.primerlist,params.fastqfile,params.samplename)

    ASSEMBLY_ch = FLYEASSEMBLY(SPLITREADS_ch)

    
}
