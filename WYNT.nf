#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include {FLYEASSEMBLY} from "./GENOTYPING.nf"



// Main workflow script for the pipeline

workflow{

    FLYEASSEMBLY()

    
}
