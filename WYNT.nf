#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include {FLYEASSEMBLY;MINISAM;HAPDUP} from "./GENOTYPING.nf"
include {SPLITTER} from "./Splitter.nf"
include {BLASTN; HLAGENOTYPER; MAKEBLASTDB;CATBLAST;BLASTNC} from "./BLAST.nf"
include {CLAIR3; VCFTOFASTA} from "./VariantCalling.nf"
include {AVA; ISONCLUST; CLUSTERSPLITTER; CLUSTERALIGNER} from "./MSA.nf"

// Main workflow script for the pipeline

workflow{

    Channel.fromPath(params.primerlist).set{primer_ch}
    Channel.fromPath(params.fastqfile).set{fastq_ch}
    
    SPLITREADS_ch = SPLITTER(primer_ch,fastq_ch,params.samplename)


    ClusterAlignments_ch = CLUSTERALIGNER(SPLITREADS_ch[0], SPLITREADS_ch[1].flatten())
    Clusters_ch = CLUSTERSPLITTER(ClusterAlignments_ch)
    //Clusters_ch[1].view()
    
    //Clusters_ch[1].flatten().view()
    Clusters_ch.transpose().view()
    //Clusters_ch.tranpose().view()
    ASSEMBLY_ch = FLYEASSEMBLY(Clusters_ch.transpose())
    //ASSEMBLY_ch = FLYEASSEMBLY(SPLITREADS_ch[0], SPLITREADS_ch[1].flatten())
    //AVA_ch = AVA(SPLITREADS_ch[0], SPLITREADS_ch[1].flatten())
    
    
    //MINISAM_ch = MINISAM(ASSEMBLY_ch)
    //CLAIR3_ch = CLAIR3(MINISAM_ch)
    //VCFTOFASTA_ch = VCFTOFASTA(CLAIR3_ch)

    //HAPDUP_ch = HAPDUP(MINISAM_ch)
    //PEPPER_ch = PEPPER(MINISAM_ch)
    BLAST_ch = BLASTN(ASSEMBLY_ch).view()
    
    BLASTcat_ch = CATBLAST(BLAST_ch)
    
    
    

    HLA_type_ch = HLAGENOTYPER(BLASTcat_ch)
}
