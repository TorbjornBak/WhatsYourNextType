#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include {FLYEASSEMBLY;FLYEASSEMBLY2;MINISAM;HAPDUP} from "./GENOTYPING.nf"
include {SPLITTER} from "./Splitter.nf"
include {BLASTN; HLAGENOTYPER; MAKEBLASTDB;CATBLAST} from "./BLAST.nf"
include {CLAIR3; VCFTOFASTA} from "./VariantCalling.nf"
include {AVA} from "./MSA.nf"
include {SPLITTER2} from "./PreAssemblyClustering.nf"

// Main workflow script for the pipeline

workflow{

    Channel.fromPath(params.primerlist).set{primer_ch}
    Channel.fromPath(params.fastqfile).set{fastq_ch}
    
    SPLITREADS_ch = SPLITTER(primer_ch,fastq_ch,params.samplename)

    ASSEMBLY_ch = FLYEASSEMBLY(SPLITREADS_ch[0], SPLITREADS_ch[1].flatten())

    SPLITTER2_ch = SPLITTER2(ASSEMBLY_ch)

    ASSEMBLY2_ch = FLYEASSEMBLY2(SPLITTER2_ch)

//    AVA_ch = AVA(SPLITREADS_ch[0], SPLITREADS_ch[1].flatten())

    
//    MINISAM_ch = MINISAM(ASSEMBLY_ch)
//    CLAIR3_ch = CLAIR3(MINISAM_ch)
    //VCFTOFASTA_ch = VCFTOFASTA(CLAIR3_ch)

    //HAPDUP_ch = HAPDUP(MINISAM_ch)
    //PEPPER_ch = PEPPER(MINISAM_ch)

    BLAST_ch = BLASTN(ASSEMBLY2_ch).groupTuple(size: 7)
    BLASTcat_ch = CATBLAST(BLAST_ch)

    
    

    HLA_type_ch = HLAGENOTYPER(BLASTcat_ch)
}
