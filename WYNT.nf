#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include {FLYEASSEMBLY;MINISAM;HAPDUP;PEPPER} from "./GENOTYPING.nf"
include {SPLITTER} from "./Splitter.nf"
include {BLASTN; HLAGENOTYPER; MAKEBLASTDB;CATBLAST} from "./BLAST.nf"


// Main workflow script for the pipeline

workflow{

    // if (params.makeblastdb == "None") {

    Channel.fromPath(params.primerlist).set{primer_ch}
    Channel.fromPath(params.fastqfile).set{fastq_ch}
    
    SPLITREADS_ch = SPLITTER(primer_ch,fastq_ch,params.samplename)

    ASSEMBLY_ch = FLYEASSEMBLY(SPLITREADS_ch[0], SPLITREADS_ch[1].flatten())
    
    //MINISAM_ch = MINISAM(ASSEMBLY_ch)
    //HAPDUP_ch = HAPDUP(MINISAM_ch)
    //PEPPER_ch = PEPPER(MINISAM_ch)
    BLAST_ch = BLASTN(ASSEMBLY_ch).groupTuple(size: 8)
    BLAST_ch.view()
    BLASTcat_ch = CATBLAST(BLAST_ch).view()

    
    
    

    HLA_type_ch = HLAGENOTYPER(BLASTcat_ch)
    // }
    
    // else {
    //     MAKEBLASTDB(params.makeblastdb,params.blastdbfastafile)
    // }
    

}
