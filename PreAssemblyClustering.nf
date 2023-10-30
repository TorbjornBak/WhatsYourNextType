process SPLITTER2{


    cpus 8
    memory '4 GB'
    

    
    input:
    tuple val(sample_name), val(allelename), path(assembly), path(splitted_reads)
    

    output:
    tuple val(sample_name), val(allelename), path("Clusters_${allelename}/*")
    
    
    
    
    script:
    """
    rm -rf -f Clusters_${allelename}
    mkdir Clusters_${allelename}
    python ${projectDir}/Clustering.py ${assembly}/assembly.fasta ${allelename} ${splitted_reads} \$PWD
    """
}


