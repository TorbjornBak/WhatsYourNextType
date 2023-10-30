#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include {FLYEASSEMBLY;MINISAM;HAPDUP;SHASTA} from "./GENOTYPING.nf"
include {SPLITTER} from "./Splitter.nf"
include {BLASTN; HLAGENOTYPER; MAKEBLASTDB;CATBLAST;BLASTNC} from "./BLAST.nf"
include {CLAIR3; VCFTOFASTA} from "./VariantCalling.nf"
include {AVA; ISONCLUST; CLUSTERSPLITTER; CLUSTERALIGNER} from "./MSA.nf"

// Main workflow script for the pipeline

workflow{

    //Channel.fromPath(params.primerlist).set{primer_ch}
    
    Channel.fromPath(params.fastqfile).set{fastq_ch}
    
    SPLITREADS_ch = SPLITTER(fastq_ch,params.samplename)
    ASSEMBLY_ch = FLYEASSEMBLY(SPLITREADS_ch.transpose())
    
    MINISAM_ch = MINISAM(ASSEMBLY_ch)

    HAPDUP_ch = HAPDUP(MINISAM_ch)
 
    BLAST_ch = BLASTN(HAPDUP_ch.transpose()).groupTuple().view()
    
    BLASTcat_ch = CATBLAST(BLAST_ch)
    
    

    HLA_type_ch = HLAGENOTYPER(BLASTcat_ch)
}
