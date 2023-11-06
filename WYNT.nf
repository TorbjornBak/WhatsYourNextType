#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include {FLYEASSEMBLY;MINISAM;HAPDUP;SHASTA;DOWNSAMPLING} from "./GENOTYPING.nf"
include {SPLITTER;FIRSTDOWNSAMPLING} from "./Splitter.nf"
include {BLASTN; HLAGENOTYPER; MAKEBLASTDB;CATBLAST;BLASTNC} from "./BLAST.nf"
include {CLAIR3; VCFTOFASTA} from "./VariantCalling.nf"
include {AVA; ISONCLUST; CLUSTERSPLITTER; CLUSTERALIGNER} from "./MSA.nf"

// Main workflow script for the pipeline

workflow{
    
    FASTQ_ch = Channel.fromPath(params.fastqfile)

    FIRSTDOWNSAMPLING_ch = FIRSTDOWNSAMPLING(FASTQ_ch)
    
    SPLITREADS_ch = SPLITTER(FIRSTDOWNSAMPLING_ch)

    DOWNSAMPLED_READS_ch = DOWNSAMPLING(SPLITREADS_ch.transpose())

    ASSEMBLY_ch = FLYEASSEMBLY(DOWNSAMPLED_READS_ch)
    
    MINISAM_ch = MINISAM(ASSEMBLY_ch)

    HAPDUP_ch = HAPDUP(MINISAM_ch)
 
    BLAST_ch = BLASTN(HAPDUP_ch.transpose()).groupTuple().view()
    
    BLASTcat_ch = CATBLAST(BLAST_ch)
    
    
    

    HLA_type_ch = HLAGENOTYPER(BLASTcat_ch)
}


