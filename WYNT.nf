#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include {ASSEMBLY;MINISAM;PHASING;UNPHASED;EXTRACTMARGIN;POLISHING} from "./Scripts/GENOTYPING.nf"
include {DEMULTIPLEXING;DOWNSAMPLING_1;DOWNSAMPLING_2} from "./Scripts/Splitter.nf"
include {BLASTN; HLAGENOTYPER; CATBLAST;} from "./Scripts/BLAST.nf"


// Main workflow script for the pipeline

workflow{
    
    FASTQ_ch = Channel.fromPath(params.fastqfile)

    FIRSTDOWNSAMPLING_ch = DOWNSAMPLING_1(FASTQ_ch)
    
    DEMULTIPLEXED_READS_ch = DEMULTIPLEXING(FIRSTDOWNSAMPLING_ch)

    DOWNSAMPLED_READS_ch = DOWNSAMPLING_2(DEMULTIPLEXED_READS_ch.transpose())

    ASSEMBLY_ch = ASSEMBLY(DOWNSAMPLED_READS_ch[0])
    //ASSEMBLY_ch[0].splitFasta(record: [id: true ]).count().view()
    //ASSEMBLY_ch[1].toInteger().view()
    //if (ASSEMBLY_ch[1].toInteger() > 0) {

    ASSEMBLYSPLT_ch = ASSEMBLY_ch.splitFasta(elem:4).groupTuple(by:[0,1,2,3])
       
    ASSEMBLYSPLT_ch.branch{
        small: it[4].size() == 1
        large: it[4].size() > 1
     }.set{result}
    

    MINISAM_ch = MINISAM(result.small)

    PHASED_ch = PHASING(MINISAM_ch)
    
    UNPHASED_ch = UNPHASED(result.large)

    BLAST_INPUT_ch = UNPHASED_ch.concat(PHASED_ch[0])
    
    //BLAST_INPUT_ch.join(DOWNSAMPLED_READS_ch[1], by: [0,1]).view()

    //POLISHED_ch = POLISHING(BLAST_INPUT_ch.join(DOWNSAMPLED_READS_ch[1], by: [0,1]).transpose())

    BLAST_ch = BLASTN(BLAST_INPUT_ch.transpose()).groupTuple()

    
    BLASTcat_ch = CATBLAST(BLAST_ch)

    MARGINLOGS_ch = PHASED_ch[1].groupTuple()



    HAPLOTYPE_DIST_ch = EXTRACTMARGIN(MARGINLOGS_ch)

    BLASTcat_ch.join(HAPLOTYPE_DIST_ch).view()

    HLA_type_ch = HLAGENOTYPER(BLASTcat_ch.join(HAPLOTYPE_DIST_ch))
}


