#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include {FLYEASSEMBLY;MINISAM;HAPDUP;PEPPER} from "./GENOTYPING.nf"
include {SPLITTER} from "./Splitter.nf"
include {BLASTN} from "./BLAST.nf"


// Main workflow script for the pipeline

workflow{

    Channel.fromPath(params.primerlist).set{primer_ch}
    Channel.fromPath(params.fastqfile).set{fastq_ch}
    
    SPLITREADS_ch = SPLITTER(primer_ch,fastq_ch,params.samplename)

    ASSEMBLY_ch = FLYEASSEMBLY(SPLITREADS_ch[0], SPLITREADS_ch[1].flatten())
    
    MINISAM_ch = MINISAM(ASSEMBLY_ch)
    //HAPDUP_ch = HAPDUP(MINISAM_ch)
    PEPPER_ch = PEPPER(MINISAM_ch)
    BLAST_ch = BLASTN(ASSEMBLY_ch)
    

}
