#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include {FASTERQDUMP;TRIM; KRAKEN; TAXREMOVE; FASTQC} from "./Trimming.nf"



// Main workflow script for the pipeline

workflow{


    
}
