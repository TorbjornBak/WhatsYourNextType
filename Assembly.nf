#!/usr/bin/env nextflow

// Script that contains nextflow processes for assembling Illumina reads, as well as for detecting the offset of 
// the reads and calculating a N50 score.

// Assembly using spades
process SPADES{
    if (params.server) {
        beforeScript 'module load spades/3.15.5'
        cpus {16 * task.attempt}
        memory { 16.GB + 32.GB * task.attempt }
        //If you change memory for spades, remember to also change the memory limit tag in the spades command below!
        errorStrategy { task.attempt < 4 ? 'retry' : 'ignore' }
        maxRetries 4
        afterScript 'module unload spades/3.15.5'
        time { 12.hour * task.attempt + (1.hours *(task.attempt - 1))}
        
    }
    else {
        conda "spades=3.15.5 conda-forge::openmp"
        cpus 8
        memory '4 GB'
        time { 1.hour * task.attempt }
    }
    
    
    publishDir "${params.outdir}/${pair_id}", mode: 'copy'

    
    input: 
    val(pair_id)
    path(r1)
    path(r2)

    val phred

    output:
    tuple val(pair_id), path ("Assembly/${params.contigs}.fasta.gz")

    
    script:
    """
    gzip -d -f ${r1}
    gzip -d -f ${r2}
    spades.py -o Assembly -1 ${r1.baseName} -2 ${r2.baseName} --meta --threads ${task.cpus} --memory ${ 16 + 32*task.attempt} --phred-offset ${phred} 
    gzip -n Assembly/${params.contigs}.fasta   
    gzip ${r1.baseName}
    gzip ${r2.baseName}
    """
}


// Offset detection using offsetdetector.py

process OFFSETDETECTOR{
    cpus 2
    time = 40.m
    memory 1.GB
    input:
    tuple(val(pair_id), path(reads))


    output:
    stdout

    script:
    
    """
    gzip -d ${reads[0]} -f
    gzip -d ${reads[1]} -f
    python3 ${projectDir}/offsetdetector.py ${reads[0].baseName} ${reads[1].baseName}
    rm ${reads[0].baseName}
    rm ${reads[1].baseName}
    """


}


// Calculating N50 score from contigs using the bbmap stats.sh script

process N50{
    if (params.server) {
        beforeScript 'module load bbmap'
        cpus 1
        time = 5.m
        memory '1 GB'
    }
    else {
        conda "agbiome::bbtools"
        cpus 1
        memory '1 GB'
    }
    publishDir "${params.outdir}/${pair_id}", mode: 'copy'



    input: 
    tuple val (pair_id), path(contigs_fasta)
    


    output:    

    script:
    """
    gzip -f -d ${contigs_fasta}
    stats.sh in=${contigs_fasta.baseName} >> ${projectDir}/${params.outdir}/${pair_id}/Assembly/assemblyStats.txt
    gzip ${contigs_fasta.baseName}
    """

}
