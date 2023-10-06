process SPLITTER2{


    cpus 8
    memory '4 GB'
    
    publishDir "${params.ouASSEMBLY_chtdir}/${samplename}", mode: 'copy'

    
    input:
    tuple val(sample_name), val(allelename), path(assembly), path(splitted_reads)
    

    output:
    tuple val(sample_name), val(allelename), path(splitted_reads), path("${allelename}.cluster1*"),path("${allelename}.cluster2*")
    
    
    
    script:
    """
    minimap2 -x ava-ont ${splitted_reads} ${assembly}/assembly.fasta > output.paf
    python ${projectDir}/Clustering.py output.paf ${allelename} ${splitted_reads}
    """
}