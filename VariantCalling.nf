#!/usr/bin/env nextflow

process CLAIR3{
    
    conda "bioconda::clair3 python=3.9.0"
    
    cpus 8
    memory '4 GB'
    time 1.hour
        
    publishDir "${params.outdir}/${sample_name}", mode: 'copy'

    
    input: 
    tuple val(sample_name), val(allelename), path(assembly), path(lr_mapping_bam), path(lr_mapping_bam_bai)
    

    output:
    tuple val(sample_name), val(allelename), path(assembly), path("${sample_name}_${allelename}_variant")
    
    script:
    """    
 
    model_name=\$CONDA_PREFIX
    model_name+="/bin/models/r941_prom_sup_g5014"
    echo \$model_name

    run_clair3.sh --bam_fn=${lr_mapping_bam} \
    --ref_fn=${assembly}/assembly.fasta --threads=${task.cpus} \
    --platform="hifi" --model_path=\$model_name \
    --output=${sample_name}_${allelename}_variant --enable_phasing \
    --include_all_ctgs --gvcf
    """
}



process VCFTOFASTA{
    
    conda "bioconda::vcflib"
    
    cpus 8
    memory '4 GB'
    time 1.hour
        
    publishDir "${params.outdir}/${sample_name}", mode: 'copy'

    
    input: 
    tuple val(sample_name), val(allelename), path(assembly), path(variantfolder)
    

    output:
    tuple val(sample_name), val(allelename), path("${sample_name}_fa*")
    
    script:
    """    
    gzip -d ${variantfolder}/merge_output.vcf.gz
    vcf2fasta --reference ${assembly} -p ${sample_name}_fa ${variantfolder}/merge_output.vcf
    """
}