#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include {FLYEASSEMBLY;MINISAM;HAPDUP;SHASTA;DOWNSAMPLING;UNPHASED} from "./GENOTYPING.nf"
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
    //ASSEMBLY_ch[0].splitFasta(record: [id: true ]).count().view()
    //ASSEMBLY_ch[1].toInteger().view()
    //if (ASSEMBLY_ch[1].toInteger() > 0) {

    ASSEMBLYSPLT_ch = ASSEMBLY_ch.splitFasta(elem:4).groupTuple(by:[0,1,2,3])
       
    ASSEMBLYSPLT_ch.branch{
        small: it[4].size() == 1
        large: it[4].size() > 1
     }.set{result}
    

    MINISAM_ch = MINISAM(result.small)

    HAPDUP_ch = HAPDUP(MINISAM_ch)
    
    UNPHASED_ch = UNPHASED(result.large)

    BLASTIN_ch = UNPHASED_ch.concat(HAPDUP_ch)    

    BLAST_ch = BLASTN(BLASTIN_ch.transpose()).groupTuple()
    
    BLASTcat_ch = CATBLAST(BLAST_ch)

    HLA_type_ch = HLAGENOTYPER(BLASTcat_ch)
}


